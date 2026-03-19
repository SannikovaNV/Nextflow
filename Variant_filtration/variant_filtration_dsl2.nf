#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ================================================================
//  V A R I A N T  F I L T R A T I O N  — DSL2  v2
// ================================================================

def helpMessage() {
    log.info """
    ================================================================
    ${workflow.manifest.name}  ~  version ${workflow.manifest.version}
    ================================================================

    Usage:
        nextflow run variant_filtration.nf [OPTIONS]

    Options:
      --vcf_files         Glob pattern for input VCF + index pairs (default: results/raw_variant_calling/*.{vcf,vcf.idx})
      --reference         Path to reference FASTA (required)
      --region_intervals  BED file for target regions — optional, WES/panels only
      --outdir            Output directory (default: results)
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// ----------------------------------------------------------------
// Default parameters
// ----------------------------------------------------------------
params.vcf_files        = "results/raw_variant_calling/*.{vcf,vcf.idx}"
params.reference        = null
params.region_intervals = "NO_FILE"
params.outdir           = "results"

if (!params.reference) {
    error "Reference FASTA is required. Use --reference /path/to/reference.fa"
}

log.info """
================================================================
V A R I A N T  F I L T R A T I O N  v2
================================================================
vcf_files        : ${params.vcf_files}
reference        : ${params.reference}
region_intervals : ${params.region_intervals}
outdir           : ${params.outdir}
================================================================
"""

// ----------------------------------------------------------------
// Processes
// ----------------------------------------------------------------

process SELECT_VARIANTS {
    tag "SelectVariants: ${sample_id}"
    label 'gatk'
    publishDir "${params.outdir}/intermediate/selected", mode: 'copy'

    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(sample_id), path(vcf_file), path(vcf_idx)

    output:
    tuple val(sample_id), path("${sample_id}.snp.vcf"), path("${sample_id}.snp.vcf.idx"),
                          path("${sample_id}.indel.vcf"), path("${sample_id}.indel.vcf.idx")

    script:
    """
    gatk SelectVariants \
        -V ${vcf_file} \
        -select-type SNP \
        -O ${sample_id}.snp.vcf

    gatk SelectVariants \
        -V ${vcf_file} \
        -select-type INDEL \
        -O ${sample_id}.indel.vcf
    """
}

process FILTER_VARIANTS {
    tag "VariantFiltration: ${sample_id}"
    label 'gatk'
    publishDir "${params.outdir}/intermediate/filtered", mode: 'copy'

    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(sample_id), path(snp_vcf), path(snp_idx),
                          path(indel_vcf), path(indel_idx)

    output:
    tuple val(sample_id), path("${sample_id}.snp.filtered.vcf"),
                          path("${sample_id}.indel.filtered.vcf")

    script:
    // GATK Best Practices hard filters
    """
    gatk VariantFiltration \
        -V ${snp_vcf} \
        -filter 'QD < 2.0'              --filter-name 'QD2' \
        -filter 'QUAL < 30.0'           --filter-name 'QUAL30' \
        -filter 'SOR > 3.0'             --filter-name 'SOR3' \
        -filter 'FS > 60.0'             --filter-name 'FS60' \
        -filter 'MQ < 40.0'             --filter-name 'MQ40' \
        -filter 'MQRankSum < -12.5'     --filter-name 'MQRankSum-12.5' \
        -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8' \
        -O ${sample_id}.snp.filtered.vcf

    gatk VariantFiltration \
        -V ${indel_vcf} \
        -filter 'QD < 2.0'               --filter-name 'QD2' \
        -filter 'QUAL < 30.0'            --filter-name 'QUAL30' \
        -filter 'FS > 200.0'             --filter-name 'FS200' \
        -filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20' \
        -O ${sample_id}.indel.filtered.vcf
    """
}

process MERGE_VCFS {
    tag "MergeVcfs: ${sample_id}"
    label 'gatk'
    publishDir "${params.outdir}/curated", mode: 'copy'

    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(sample_id), path(snp_filtered), path(indel_filtered)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.curated.vcf"),
                          path("${sample_id}.curated.vcf.idx")

    script:
    // MergeVcfs handles sorting internally — no need for SortVcf beforehand
    """
    gatk MergeVcfs \
        -I ${snp_filtered} \
        -I ${indel_filtered} \
        -R ${reference} \
        -O ${sample_id}.curated.vcf

    gatk IndexFeatureFile \
        --input ${sample_id}.curated.vcf
    """
}

// ----------------------------------------------------------------
// Workflow
// ----------------------------------------------------------------

workflow {

    // Reference
    ch_reference = Channel.fromPath(params.reference, checkIfExists: true)

    // Input VCFs — expects pairs: sample.vcf + sample.vcf.idx
    ch_vcf = Channel
        .fromFilePairs(params.vcf_files, size: 2, checkIfExists: true)
        .map { sample_id, files ->
            def vcf = files.find { it.name.endsWith('.vcf') }
            def idx = files.find { it.name.endsWith('.idx') }
            tuple(sample_id, vcf, idx)
        }

    // Pipeline
    ch_selected = SELECT_VARIANTS(ch_vcf)
    ch_filtered = FILTER_VARIANTS(ch_selected)
    MERGE_VCFS(ch_filtered, ch_reference)
}
