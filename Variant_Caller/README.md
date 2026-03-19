# Variant Caller Pipeline

A Nextflow DSL2 pipeline for germline variant calling from paired-end FASTQ reads. Supports two aligners (BWA, Bowtie2) and three variant callers (GATK HaplotypeCaller, FreeBayes, VarScan), configurable at runtime. Designed for both local execution and HPC clusters via SLURM, with Docker and Singularity support.

---

## Pipeline overview

```
Paired-end FASTQ reads
        │
        ▼
Reference indexing     ── BWA index or Bowtie2 index
        │                  GATK sequence dictionary
        │                  SAMtools faidx
        ▼
Alignment              ── BWA-MEM or Bowtie2
        │
        ▼
Mark duplicates        ── GATK MarkDuplicates (optional removal)
        │
        ▼
BAM indexing           ── SAMtools index
        │
        ▼
Variant calling        ── GATK HaplotypeCaller
                          FreeBayes
                          VarScan (SNPs + indels separately)
```

---

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 23.0.0
- Docker or Singularity

No other local installations required — all tools run inside containers.

---

## Usage

**Run locally with Docker:**
```bash
nextflow run variant_caller.nf \
    -profile docker \
    --reads 'data/*R{1,2}*.fastq.gz' \
    --genome /path/to/reference.fa \
    --outdir results
```

**Run on HPC with SLURM + Singularity:**
```bash
nextflow run variant_caller.nf \
    -profile slurm,singularity \
    --reads 'data/*R{1,2}*.fastq.gz' \
    --genome /path/to/reference.fa \
    --outdir results
```

**Use FreeBayes instead of GATK:**
```bash
nextflow run variant_caller.nf \
    -profile docker \
    --reads 'data/*R{1,2}*.fastq.gz' \
    --genome /path/to/reference.fa \
    --vc freebayes \
    --outdir results
```

**WES or panel — restrict to target regions:**
```bash
nextflow run variant_caller.nf \
    -profile docker \
    --reads 'data/*R{1,2}*.fastq.gz' \
    --genome /path/to/reference.fa \
    --region_intervals /path/to/targets.bed \
    --outdir results
```

**Alignment only — skip variant calling:**
```bash
nextflow run variant_caller.nf \
    -profile docker \
    --reads 'data/*R{1,2}*.fastq.gz' \
    --genome /path/to/reference.fa \
    --skip_variant_calling true \
    --outdir results
```

---

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--reads` | Glob pattern for paired-end FASTQ files | `data/*R{1,2}*.fastq.gz` |
| `--genome` | Path to reference FASTA | required |
| `--outdir` | Output directory | `results` |
| `--region_intervals` | BED file for target regions (WES/panels) | off |
| `--aln` | Aligner: `bwa` or `bowtie2` | `bwa` |
| `--vc` | Variant caller: `gatk`, `freebayes`, `varscan` | `gatk` |
| `--remove_duplicates` | Remove marked duplicates | `false` |
| `--skip_variant_calling` | Align only, skip variant calling | `false` |
| `--min_alt_fraction` | Minimum allele frequency (FreeBayes, VarScan) | `0.2` |
| `--ploidy` | Sample ploidy | `2` |
| `--threads` | Threads per process | `4` |
| `--platform` | Sequencing platform for read group | `ILLUMINA` |

---

## Output

```
results/
├── reference/         # Index files and sequence dictionary
├── alignment/         # Sorted, deduplicated and indexed BAM files
└── vcf/               # Variant calls per sample
    ├── *.gatk.vcf         # GATK HaplotypeCaller output
    ├── *.freebayes.vcf    # FreeBayes output
    ├── *.varscan.snp.vcf  # VarScan SNP calls
    └── *.varscan.indel.vcf # VarScan indel calls
```

---

## Tools and containers

| Tool | Version | Container |
|------|---------|-----------|
| BWA | 0.7.17 | `biocontainers/bwa:v0.7.17` |
| Bowtie2 | 2.3.5.1 | `biocontainers/bowtie2:v2.3.5.1` |
| SAMtools | 1.9 | `biocontainers/samtools:v1.9` |
| GATK | 4.1.4.1 | `broadinstitute/gatk:4.1.4.1` |
| FreeBayes | 1.3.1 | `biocontainers/freebayes:v1.3.1` |
| VarScan | 2.4.4 | `biocontainers/varscan:v2.4.4` |

---

## Notes

- Reference indexing runs once per execution. If the index files already exist in the output directory, point `--genome` to the pre-indexed reference to skip this step.
- VarScan produces two separate VCF files per sample — one for SNPs and one for indels. These can be merged downstream with `bcftools merge` if needed.
- The `--region_intervals` BED file is optional. When not provided the pipeline runs in WGS mode across the entire genome.
- Resume a failed run from the last successful step with `-resume`.

---

## Author

Natalia Sannikova
