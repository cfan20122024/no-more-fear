# RNA-seq Quantification Workflow

**Workflow Name:** `rna-seq-quant`  
**Workflow Language:** Nextflow DSL2  
**Latest Release:** v3.0.4 (as of 2025-11-12)  
**Repository:** UNC LBG Workflows – RNA-seq Quantification  
**Workflow URL:** [UNC LBG Workflows – RNA-seq Quantification](https://sc.unc.edu/lbg/workflows/nextflow/rna-seq-quant/-/tags)

## Workflow Summary
This Nextflow DSL2 workflow performs alignment, quantification, and quality control of RNA-seq data generated from Illumina platforms. It supports both paired-end and single-end FASTQ inputs and produces:

- Gene- and transcript-level expression tables  
- Normalized expression matrices (TPM and upper-quartile)  
- Comprehensive quality control reports

Input files are provided in a **tab-separated values (TSV) file** that lists sample IDs and paths to FASTQ files. Multiple FASTQ files per sample can be merged automatically.

## Key Features
- Automated end-to-end processing from raw FASTQs to final normalized expression outputs  
- Integration of multiple QC tools: FastQC, fastp, Picard, samtools, and MultiQC  
- Generation of TPM and upper-quartile normalized (quantile75eq1k) expression matrices  
- Compatibility with standard human and mouse reference transcriptomes (e.g., Ensembl GRCh38, GRCh38/Gencode, mm10)  
- Modular design using Nextflow DSL2 sub-workflows and containerization for reproducibility  
- Transparent versioning and reproducibility: workflow version is determined by the release/tag

## Workflow Versioning
- Latest workflow version as of November 12, 2025: **v3.0.4**  
- Version is explicitly set via Nextflow `-r` parameter  
- Users can verify tool versions automatically logged during the run (`rnaseq_mqc_versions.yml`) or in MultiQC reports  
- Full list of available workflow versions: [https://sc.unc.edu/lbg/workflows/nextflow/rna-seq-quant/-/tags](https://sc.unc.edu/lbg/workflows/nextflow/rna-seq-quant/-/tags)

## Workflow Outline

1. **Input Preparation**
   - Accepts a **tab-separated sample sheet (TSV)** listing sample IDs and paths to FASTQ files.  
   - Multiple FASTQ files per sample (e.g., from different lanes) are automatically merged.

2. **Read Preprocessing**
   - Trim adapters and low-quality bases using **fastp v0.24.0**.  
   - Generates per-sample **pre- and post-trimming statistics**.

3. **Alignment**
   - Map reads to the reference genome using **STAR v2.7.11b**.

4. **Quantification**
   - Quantify gene and transcript abundances using **Salmon v1.10.3** (alignment-based mode).  
   - Outputs **raw counts** and **TPM values**. Most downstream analyses use either TPM or Salmon counts.

5. **Normalization**
   - Compute **upper-quartile normalized expression values** (`quantile75eq1k`) using a custom R script (**quantile75to1k.R**).  
     - Quantile (0.75) scaled to 1000.  
   - Both **TPM** and **quantile75eq1k matrices** are provided for downstream analysis.

6. **Quality Control (QC)**
   - **FastQC v0.12.1**: read-level quality assessment.  
   - **fastp v0.24.0**: trimming statistics.  
   - **Picard v2.22.4**: `CollectRnaSeqMetrics` for alignment metrics.  
   - **samtools v1.10**: `flagstat` for mapping statistics.  
   - Custom QC metrics: maximum read length, total reads, mapped/unmapped counts.

7. **Aggregation and Reporting**
   - Merge expression matrices across all samples.  
   - Summarize QC metrics using **MultiQC v1.27.1**.  
   - Generate a **custom HTML QC summary** (`RNA_QC_report.html`) across all samples.

8. **Annotation**
   - Genes are annotated by **Ensembl Gene ID**, **Entrez Gene ID**, and **Gene Name**.  
   - Note: Gene name outputs may contain duplicates if multiple Ensembl IDs share the same symbol.

## Tools and Versions

| Step             | Tool                | Version   | Purpose / Notes                               |
|-----------------|-------------------|----------|-----------------------------------------------|
| Read trimming    | fastp              | 0.24.0   | Adapter and quality trimming, QC stats        |
| Alignment        | STAR               | 2.7.11b  | Spliced read alignment                        |
| Quantification   | Salmon             | 1.10.3   | Transcript-level quantification               |
| QC               | FastQC             | 0.12.1   | Read quality summaries                         |
| QC               | Picard             | 2.22.4   | RNA-seq metrics collection (CollectRnaSeqMetrics) |
| QC               | samtools           | 1.10     | Alignment statistics                           |
| QC summary       | MultiQC            | 1.27.1   | Aggregated QC reporting                        |
| Normalization    | quantile75to1k.R   | custom   | Upper-quartile scaling normalization           |

### Version Verification
- Tool versions are automatically logged during each run in the `rnaseq_mqc_versions.yml` file.  
- Verified versions correspond to workflow **v3.0.4** (as tested November 2025).  
- Alternatively, view the “Software Version” section in the MultiQC report.

## Input Files

The RNA-seq quant workflow requires a **tab-separated values (TSV) file** listing the input FASTQ files for each sample. The format depends on whether your data is **paired-end** or **single-end**, and whether you need to merge multiple FASTQ files per sample.

---

### 1. Running in jUNCtion

When running the workflow in jUNCtion, you must provide a TSV file describing your samples. The workflow reads this file to locate your FASTQ files.

**Requirements for the TSV file:**
- Must be tab-separated.
- Must include a header row.
- Columns:
  - `sample`: Unique sample identifier (used in output filenames)
  - `fq1`: Path to the first FASTQ file (R1)
  - `fq2`: Path to the second FASTQ file (R2) for paired-end data. Omit for single-end data.

**Notes for jUNCtion:**
- You can provide multiple rows for a single sample if you need to merge multiple FASTQ files (e.g., multiple lanes).
- The workflow automatically combines FASTQ files for the same sample.

---

### 2. Running on the Command Line

The TSV file format is identical on the command line. The workflow reads the TSV file using the `--input` parameter:

```bash
nextflow run rna_seq_quant.nf --input path/to/samples.tsv
```

## Output Files

### 1. Salmon Quantification Files

The Salmon workflow produces gene- and isoform-level expression estimates. Files are available in both **raw** and **normalized** forms, with identifiers as Ensembl Gene ID, Entrez Gene ID, or Gene Symbol.

| File | Description | Data Level | Identifier | Generated By |
|------|------------|-----------|-----------|--------------|
| salmon_ensemblgid.matrix.tsv | Raw Salmon expression estimates per gene | Gene | Ensembl Gene ID | `collect_column.pl` aggregates `quant.sf` from each sample |
| salmon_ensemblgid_quantile75eq1k.matrix.tsv | Upper-quartile normalized Salmon expression | Gene | Ensembl Gene ID | Custom normalization script |
| salmon_entrez_gene.matrix.tsv | Raw Salmon expression | Gene | Entrez Gene ID | `collect_column.pl` |
| salmon_entrez_gene_quantile75eq1k.matrix.tsv | Upper-quartile normalized Salmon expression | Gene | Entrez Gene ID | Custom normalization script |
| salmon_gene.matrix.tsv | Raw Salmon expression | Gene | Gene Symbol | `collect_column.pl` |
| salmon_gene_quantile75eq1k.matrix.tsv | Upper-quartile normalized Salmon expression | Gene | Gene Symbol | Custom normalization script |
| salmon_isoform.matrix.tsv | Salmon expression at isoform (transcript) level | Transcript | Ensembl Transcript ID | `collect_column.pl` |
| salmon_tpm_ensemblgid.matrix.tsv | TPM (Transcripts Per Million) estimates | Gene | Ensembl Gene ID | `collect_column.pl` |
| salmon_tpm_entrez_gene.matrix.tsv | TPM estimates | Gene | Entrez Gene ID | `collect_column.pl` |
| salmon_tpm_gene.matrix.tsv | TPM estimates | Gene | Gene Symbol | `collect_column.pl` |
| salmon_tpm_isoform.matrix.tsv | TPM estimates at isoform level | Transcript | Ensembl Transcript ID | `collect_column.pl` |

**Notes:**
- Raw Salmon outputs for each sample are stored in the `salmon.quant` folder.
- The `collect_column.pl` Perl script extracts the `quant.sf` column from each sample’s Salmon output and collates them into cohort-level matrices.

---

### 2. Per-Sample Output Files (--mode sample)

These files are generated individually for each sample.

| File Type | Description | Generated By |
|-----------|------------|--------------|
| `*.sort.bam` | Aligned and coordinate-sorted BAM file (if BAM publishing is enabled) | STAR aligner |
| `*.fastqc.html` | FastQC quality control report | FastQC |
| `*.fastp.html` | Fastp quality control report | Fastp |
| `*.samtools.txt` | Alignment statistics (flagstat) | Samtools |
| `*.rnaqc.txt` | RNA-seq QC metrics | Picard CollectRnaSeqMetrics |
| `*salmon_*.txt` | Per-sample Salmon quantification files | `collect_column.pl` collates `quant.sf` |
| `salmon.quant/` | Raw Salmon quantification folder for each sample | Salmon |

---

### 3. Cohort Summary Files (--mode summary)

Aggregated files summarizing all samples in the cohort.

#### 3.1 Salmon Expression Matrices

- `salmon_ensemblgid.matrix.tsv`: Raw gene-level counts with Ensembl IDs
- `salmon_ensemblgid_quantile75eq1k.matrix.tsv`: Upper-quartile normalized gene counts with Ensembl IDs
- `salmon_entrez_gene.matrix.tsv`: Raw gene-level counts with Entrez Gene IDs
- `salmon_entrez_gene_quantile75eq1k.matrix.tsv`: Upper-quartile normalized gene counts with Entrez Gene IDs
- `salmon_gene.matrix.tsv`: Raw gene-level counts with Gene Symbols
- `salmon_gene_quantile75eq1k.matrix.tsv`: Upper-quartile normalized gene counts with Gene Symbols
- `salmon_isoform.matrix.tsv`: Isoform-level raw counts
- `salmon_tpm_ensemblgid.matrix.tsv`: TPM estimates with Ensembl Gene IDs
- `salmon_tpm_entrez_gene.matrix.tsv`: TPM estimates with Entrez Gene IDs
- `salmon_tpm_gene.matrix.tsv`: TPM estimates with Gene Symbols
- `salmon_tpm_isoform.matrix.tsv`: TPM isoform-level estimates

All Salmon matrices are generated using `collect_column.pl` to collate individual sample `quant.sf` files, with additional normalization steps applied where indicated.

#### 3.2 MultiQC Report

- `multiqc_report.html`: Comprehensive QC summary across all samples  
  Combines metrics from FastQC, Fastp, Samtools, Picard, and Salmon  
  Generated by MultiQC

#### 3.3 Custom RNA QC Report

- `RNA_QC_report.html`: Detailed QC report generated by `RNA_QC_report.Rmd`  
  Includes:
  - Number of samples in the dataset
  - Number of genes with counts >0 and >4
  - Relative Log Expression (RLE) boxplots for normalization
  - Principal Component Analysis (PCA) plots (PC1–PC10)
  - Picard QC plots: total reads, aligned reads, mRNA bases, % ribosomal bases
  - Boxplots of aligned counts
  - Upper-quartile normalization factor plots
  - GAPDH and top-expressed gene expression plots
  - Unsupervised heatmap
  - List of potential outliers (informational; dependent on the dataset)

#### 3.4 Filtered Counts

- `filteredCounts.rds`: Filtered Salmon count matrix (genes with max count >10) used as input for PCA and variance-stabilized transformation (vst)  
  - Format: `.rds`  
  - Generated by internal R script during QC pipeline


## Running at UNC

The RNA-seq quant workflow can be executed either via the **jUNCtion workflow manager** or locally on the **LBG cluster** using Nextflow.

- **Nextflow Documentation:** For local runs on the LBG cluster, refer to [Nextflow at UNC/LBG](https://sc.unc.edu/lbg/workflows/nextflow/lbg-nextflow-docs).  
- **Workflow Execution:** When running locally, the workflow should be executed **twice**:  
  1. **Sample mode (`--mode sample`)** – processes each sample individually.  
  2. **Summary mode (`--mode summary`)** – aggregates the per-sample outputs into cohort-level matrices and QC reports. Ensure the summary run is performed **only after** the sample mode run completes.

---

### Example: Local LBG Slurm Run (hg38, paired-end)

### Example: Local LBG Slurm Run (hg38, paired-end)

The following commands demonstrate running the RNA-seq quant workflow in **sample mode** for per-sample processing and then in **summary mode** to aggregate cohort-level results. Copy and paste the block below to execute on the LBG cluster.

```bash
# 1. Per-sample processing
nextflow run https://sc.unc.edu/lbg/workflows/nextflow/rna-seq-quant \
    -with-report report.html \
    -with-trace trace.txt \
    --gencode_to_entrez_mapping /path/to/empty.txt \
    -w working \
    -profile lbg_slurm,lbg_hg38_gencode_v36 \
    -r rna-seq-quant-3.0.4 \
    --input inventory.tsv \
    --mode sample \
    --output_dir /path/to/results \
> rna_seq_quant.log 2>&1

# 2. Cohort summary run
nextflow run https://sc.unc.edu/lbg/workflows/nextflow/rna-seq-quant \
    -with-report report.summary.html \
    -with-trace trace.summary.txt \
    --gencode_to_entrez_mapping /path/to/empty.txt \
    -w working \
    -profile lbg_slurm,lbg_hg38_gencode_v36 \
    -r rna-seq-quant-3.0.4 \
    --input /path/to/results/output_inventory.tsv \
    --mode summary \
    --output_dir /path/to/results \
> rna_seq_quant_summary.log 2>&1
```

**Key Notes:**  
- The `-profile` argument specifies the Slurm execution environment and genome annotation (hg38/Gencode v36 in this example).  
- The `-with-report` and `-with-trace` flags generate HTML reports and workflow traces for monitoring performance.  
- `working` is the workflow working directory; you can specify any path accessible on your cluster.  
- Logs are captured to `rna_seq_quant.log` and `rna_seq_quant_summary.log` for debugging and record-keeping.  

## Issues and Updates

**DSL2 Parameter Error**  
- **Problem:** The `-dsl2` flag causes an "unknown option" error.  
- **Solution:** Omit the `-dsl2` parameter when running `rna-seq-quant-3.0.4`.

**Singularity "--no-home" Error**  
- **Problem:** Workflow fails during `merge_fastqs_paired` with error: `Unknown option: --no-home`  
- **Cause:** Bug when `/dev/null` is passed as a parameter value.  
- **Fix:** Use an empty file instead of omitting the parameter: `--gencode_to_entrez_mapping /path/to/empty.txt`  
- **Note:** This is a temporary workaround until developers fix Nextflow/Singularity parameter handling.

**Home Directory Access Issue**  
- **Problem:** Singularity containers cannot access files in the default home directory: `/home/username/.nextflow/`, which isn't mounted by default.  
- **Symptoms:** Failures accessing `multiqc_config.yaml` and other assets.  
- **Temporary Fix:** `NXF_SINGULARITY_HOME_MOUNT=true nextflow run [workflow]`  
- **Permanent Solution:** `export NXF_HOME=/datastore/nextgenout5/share/labs/bioinformatics/username/.nextflow`  
- **Reason:** Default `$HOME/.nextflow` isn't mounted into Singularity containers. Using `/datastore` ensures container accessibility.  
- **Note:** Replace `username` with your actual username. This solution has been tested and resolved the issue.

# Single-End Data

For single-end data:

- **Input TSV**: Omit the `fq2` column
- **Command**: Add `--single` option

```bash
nextflow run [workflow] --single
```

# Available Profiles

## Human References
- `lbg_hg38_gencode_v22` - hg38 with GENCODE v22 annotations
- `lbg_hg38_gencode_v36` - hg38 with GENCODE v36 annotations
- `lbg_hg38_gencode_v47` - hg38 with GENCODE v47 annotations

## Mouse References
- `lbg_mm10_ensembl` - mm10 with Ensembl v84 annotations
- `lbg_GRCm38_p6_vM25` - GRCm38 with GENCODE vM25/Ensembl v100 annotations

## Rat Reference
- `lbg_rn7` - Rat rn7 (mRatBN7.2 genome assembly)
