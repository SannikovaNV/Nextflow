#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ================================================================
//  V A R I A N T  C A L L E R  — DSL2  v2
// ================================================================

def helpMessage() {
    log.info """
    ================================================================
    V A R I A N T  C A L L E R
    ================================================================

    Usage:
        nextflow run variant_caller.nf [OPTIONS]

    Options:
      --indir        DIR         Input directory with FASTQ files (default: data)
      --outdir       DIR         Output directory (default: results)
      --genome       FILE        Path to reference FASTA (default: GRCh38 from iGenomes)
      --adapters     FILE        Adapter file for Trimmomatic
      --region_intervals FILE    BED file for target regions (WES/panels)

      --paired       true|false  Paired-end or single-end mode (default: true)
      --reads        GLOB        Glob pattern for paired reads
                                 (default: data/*R{1,2}*.fastq.gz)

      --aln          bwa|bowtie2 Aligner (default: bwa)
      --vc           gatk|freebayes|varscan  Variant caller (default: gatk)

      --skip_variant_calling true|false  Only align, skip variant calling (default: false)
      --remove_duplicates    true|false  Remove marked duplicates (default: false)
      --min_alt_fraction     NUM         FreeBayes min allele frequency (default: 0.2)
      --ploidy               INT         Sample ploidy (default: 2)
      --threads              INT         Threads per process (default: 4)
      --platform             STR         Sequencing platform for read group (default: ILLUMINA)
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// ----------------------------------------------------------------
// Default parameters
// ----------------------------------------------------------------
params.indir               = "data"
params.outdir              = "results"
params.genome              = "GRCh38"
params.adapters            = "${baseDir}/adapters/TruSeq3-PE.fa"
params.region_intervals    = "NO_FILE"
params.paired              = true
params.reads               = "${params.indir}/*R{1,2}*.fastq.gz"
params.aln                 = "bwa"
params.vc                  = "gatk"
params.skip_variant_calling = false
params.remove_duplicates   = false
params.min_alt_fraction    = 0.2
params.ploidy              = 2
params.threads             = 4
params.platform            = "ILLUMINA"

log.info """
================================================================
V A R I A N T  C A L L E R  v2
================================================================
genome               : ${params.genome}
reads                : ${params.reads}
region_intervals     : ${params.region_intervals}
paired               : ${params.paired}
aligner              : ${params.aln}
variant_caller       : ${params.vc}
remove_duplicates    : ${params.remove_duplicates}
skip_variant_calling : ${params.skip_variant_calling}
ploidy               : ${params.ploidy}
threads              : ${params.threads}
outdir               : ${params.outdir}
================================================================
"""

// ----------------------------------------------------------------
// Modules
// ----------------------------------------------------------------

process BWA_INDEX {
    tag "BWA index: ${reference.name}"
    label 'bwa'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path reference

    output:
    tuple path(reference), path("${reference}.*")

    script:
    """
    bwa index ${reference}
    """
}

process BOWTIE2_INDEX {
    tag "Bowtie2 index: ${reference.name}"
    label 'bowtie2'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path reference

    output:
    tuple path(reference), path("*.bt2")

    script:
    def prefix = reference.baseName
    """
    bowtie2-build --threads ${params.threads} ${reference} ${prefix}
    """
}

process GATK_DICT {
    tag "GATK sequence dictionary: ${reference.name}"
    label 'gatk'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path reference

    output:
    path "*.dict"

    script:
    """
    gatk CreateSequenceDictionary -R ${reference}
    """
}

process SAMTOOLS_FAIDX {
    tag "SAMtools faidx: ${reference.name}"
    label 'samtools'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path reference

    output:
    path "*.fai"

    script:
    """
    samtools faidx ${reference}
    """
}

process ALIGN_BWA {
    tag "BWA-MEM: ${sample_id}"
    label 'bwa'
    publishDir "${params.outdir}/alignment", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    tuple path(reference), path(index_files)

    output:
    tuple val(sample_id), path("${sample_id}.bwa.bam")

    script:
    def rg = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:${params.platform}\\tLB:${sample_id}"
    """
    bwa mem \
        -t ${params.threads} \
        -R "${rg}" \
        ${reference} \
        ${reads[0]} ${reads[1]} \
    | samtools view -bS \
    | samtools sort -@ ${params.threads} -o ${sample_id}.bwa.bam
    """
}

process ALIGN_BOWTIE2 {
    tag "Bowtie2: ${sample_id}"
    label 'bowtie2'
    publishDir "${params.outdir}/alignment", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    tuple path(reference), path(index_files)

    output:
    tuple val(sample_id), path("${sample_id}.bowtie2.bam")

    script:
    def prefix = reference.baseName
    """
    bowtie2 \
        -p ${params.threads} \
        -x ${prefix} \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
    | samtools view -bS \
    | samtools sort -@ ${params.threads} -o ${sample_id}.bowtie2.bam
    """
}

process MARK_DUPLICATES {
    tag "MarkDuplicates: ${sample_id}"
    label 'gatk'
    publishDir "${params.outdir}/alignment", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}.markdup.bam"), path("${sample_id}.markdup.metrics")

    script:
    def remove = params.remove_duplicates ? "--REMOVE_DUPLICATES true" : ""
    """
    gatk MarkDuplicates \
        -I ${bam_file} \
        -O ${sample_id}.markdup.bam \
        -M ${sample_id}.markdup.metrics \
        ${remove}
    """
}

process BAM_INDEX {
    tag "BAM index: ${sample_id}"
    label 'samtools'
    publishDir "${params.outdir}/alignment", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file), path(metrics)

    output:
    tuple val(sample_id), path(bam_file), path("${bam_file}.bai")

    script:
    """
    samtools index -@ ${params.threads} ${bam_file}
    """
}

process VARIANT_CALLING_GATK {
    tag "GATK HaplotypeCaller: ${sample_id}"
    label 'gatk'
    publishDir "${params.outdir}/vcf", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bai_file)
    path reference
    path fai
    path dict
    path intervals

    output:
    tuple val(sample_id), path("${sample_id}.gatk.vcf"), path("${sample_id}.gatk.vcf.idx")

    script:
    def region = intervals.name != 'NO_FILE' ? "-L ${intervals}" : ""
    """
    gatk HaplotypeCaller \
        --native-pair-hmm-threads ${params.threads} \
        -R ${reference} \
        -I ${bam_file} \
        -O ${sample_id}.gatk.vcf \
        ${region}
    """
}

process VARIANT_CALLING_FREEBAYES {
    tag "FreeBayes: ${sample_id}"
    label 'freebayes'
    publishDir "${params.outdir}/vcf", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bai_file)
    path reference
    path fai
    path intervals

    output:
    tuple val(sample_id), path("${sample_id}.freebayes.vcf")

    script:
    def region = intervals.name != 'NO_FILE' ? "-t ${intervals}" : ""
    """
    freebayes \
        -f ${reference} \
        -p ${params.ploidy} \
        --min-alternate-fraction ${params.min_alt_fraction} \
        ${region} \
        ${bam_file} \
    > ${sample_id}.freebayes.vcf
    """
}

process VARIANT_CALLING_VARSCAN {
    tag "VarScan: ${sample_id}"
    label 'varscan'
    publishDir "${params.outdir}/vcf", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bai_file)
    path reference
    path fai
    path intervals

    output:
    tuple val(sample_id), path("${sample_id}.varscan.snp.vcf"), path("${sample_id}.varscan.indel.vcf")

    script:
    def region = intervals.name != 'NO_FILE' ? "-l ${intervals}" : ""
    """
    samtools mpileup \
        -f ${reference} \
        ${region} \
        ${bam_file} \
    | varscan mpileup2cns \
        --variants \
        --output-vcf 1 \
        --min-var-freq ${params.min_alt_fraction} \
        --p-value 0.05 \
        --strand-filter 1 \
    > ${sample_id}.varscan.snp.vcf

    samtools mpileup \
        -f ${reference} \
        ${region} \
        ${bam_file} \
    | varscan mpileup2indel \
        --output-vcf 1 \
        --min-var-freq ${params.min_alt_fraction} \
    > ${sample_id}.varscan.indel.vcf
    """
}

// ----------------------------------------------------------------
// Workflow
// ----------------------------------------------------------------

workflow {

    // Reference
    ch_reference = Channel.fromPath(params.genome, checkIfExists: true)

    // Intervals (optional)
    ch_intervals = params.region_intervals != 'NO_FILE'
        ? Channel.fromPath(params.region_intervals, checkIfExists: true)
        : Channel.fromPath("NO_FILE", checkIfExists: false).ifEmpty { file("NO_FILE") }

    // Reads
    ch_reads = Channel
        .fromFilePairs(params.reads, size: 2, checkIfExists: true)

    // Reference preparation
    ch_fai  = SAMTOOLS_FAIDX(ch_reference)
    ch_dict = GATK_DICT(ch_reference)

    // Alignment
    if (params.aln == 'bwa') {
        ch_index   = BWA_INDEX(ch_reference)
        ch_aligned = ALIGN_BWA(ch_reads, ch_index)
    } else if (params.aln == 'bowtie2') {
        ch_index   = BOWTIE2_INDEX(ch_reference)
        ch_aligned = ALIGN_BOWTIE2(ch_reads, ch_index)
    } else {
        error "Unknown aligner: ${params.aln}. Choose bwa or bowtie2."
    }

    // Mark duplicates and index
    ch_marked  = MARK_DUPLICATES(ch_aligned)
    ch_indexed = BAM_INDEX(ch_marked)

    // Variant calling
    if (!params.skip_variant_calling) {
        if (params.vc == 'gatk') {
            VARIANT_CALLING_GATK(
                ch_indexed,
                ch_reference,
                ch_fai,
                ch_dict,
                ch_intervals
            )
        } else if (params.vc == 'freebayes') {
            VARIANT_CALLING_FREEBAYES(
                ch_indexed,
                ch_reference,
                ch_fai,
                ch_intervals
            )
        } else if (params.vc == 'varscan') {
            VARIANT_CALLING_VARSCAN(
                ch_indexed,
                ch_reference,
                ch_fai,
                ch_intervals
            )
        } else {
            error "Unknown variant caller: ${params.vc}. Choose gatk, freebayes, or varscan."
        }
    }
}
