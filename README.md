# epicypher-dafseq-pipeline

![alt text](https://www.epicypher.com/wp-content/uploads/2024/03/nav_logo-min.png)


# Deaminase Assisted Fiber-seq (DAF-seq) Nextflow Pipeline
<h3>Disclaimer</h3>
While this repository is maintained by EpiCypher, Inc., we make no promises to troubleshoot or offer technical support. Please report any bugs encountered, however we make no promises to fix them.

---

<h3>Summary</h3>

- Purpose: This is a Nextflow based analysis pipeline focused on analyzing DAF-seq data from an ONT sequencing instrument. This pipeline will carry unaligned bams through alignment, labelling, a second alignment, and ddda to 6mA conversion. Optionally, splitting reads by strand and bigwig creation are supported. 
- Entry point: `main.nf` (DSL2 workflow).


---

<h3>Requirements</h3>

- Containers:
  - `-profile local` uses Docker.
  - `-profile slurm` uses Singularity/Apptainer. (Can also be used to support Docker)

---

<h3>Inputs</h3>

<h4>Required</h4>

- `--sample_sheet`: TSV with header and columns:
  `samp_name <TAB> bam_path <TAB> ref_name <TAB> target_name`.
  - Multiple rows can share the same `samp_name` (those BAMs will be aligned individually then merged, useful for technical sequencing replicates).\

- `--ref_sheet`: TSV with header and columns: `ref_name <TAB> fasta_path <TAB> fasta_index`.
  - Can be created using `prepare_references.sh`
  - The `ref_name` values must match between the sample sheet and reference sheet.

- `--target_sheet`: TSV with header column: `chrom <TAB> start <TAB> end <TAB> target_name`:
  - `target_name` must match the name of the target in the sample sheet

- `-profile`: Nextflow explicit parameter, determines method for job execution. Use one of the following:
  - local: executes job using local resources and Docker as container method
  - aws_env: executes job using Amazon Web Service Batch compute environment
  - slurm: executes job locally using slurm as job scheduler. Singularity or Docker can be used as container method.


Example: `inputs/sample_sheet.tsv`

| samp_name | bam_path | ref_name | target_name |
|-----------|----------|----------|----------|
| sampleA | inputs/bams/a1/sampleA.bam | hg38 | ColoReg2 |
| sampleA | inputs/bams/a2/sampleA.bam | hg38 | ColoReg2 |
| sampleB | inputs/bams/b1/sampleB.bam | chm13 | NAPA |


Example `inputs/reference_sheet.tsv`:

| ref_name | ref_fasta | ref_index |
|----------|-----------|-----------|
| hg38 | /refs/hg38/hg38.fasta | /refs/hg38/hg38.fasta.fai |
| chm13 | /refs/chm13/chm13.fasta | /refs/chm13/chm13.fasta.fai |

Example `inputs/target_sheet.tsv`
| chrom | start | end | target_name |
|----------|-----------|-----------|-----------|
| chr17 | 19446516 | 19460000 | ColoReg2 |
| chr19 | 47487634 | 47515258 | NAPA |

<h4>Optional Parameters</h4>

- `--outdir` (default `${workflow.launchDir}/results`) — top-level output
  directory.
- `--confidence_ml_val` (default `250`) — ML threshold to use for both 6mA and 5mC for `ft add-nucleosomes`
  and pileups.
- `--minimum_msp_dist` (default `10`) — MSP length filter used for pileups
  (`ftx "len(msp) > ..."`).

<h4>Optional steps</h4>

- `--split_reads` (default `false`) - split reads into +/- strand based on G>A or C>T conversions
- `--create_bigwigs` (default `false`) — create pileup TSVs and BigWigs.
- `--debug` (default `false`) — prints helpful channel `view()` messages.

---

<h3>High-level workflow steps</h3>

Short descriptions of each Nextflow job step:
- Read the sample sheet and the reference sheet and join by `ref_name`.
- `align_bams`: align each input BAM with `minimap2`.
- label reads using python script that will find all transitions and label read as either G>A to C>T, then convert the majority transitions to R ambiguous IUPAC code
- Align ambiguous and labelled reads again
- Convert transitions to 6mA calls for downstream tools
- Call nucleosomes using fibertools `add-nucleosomes`
- Split reads into on/off target
- Calculate % Conversion stats for on/off target reads
- If using `--split_reads`: split reads into +/- strand based on G>A or C>T conversion rate
- If using `--create_bigwigs`: creates methylation (6mA/5mC) and nucleosome pileups and converts them to BigWigs.

---

<h3>Published output locations</h3>

- `${outdir}/1_Aligned-NotLabelled/` — Reads that have undergone the initial minimap2 alignment.
- `${outdir}/2_Aligned-Labelled/` — Reads that have undergone C>T or G>A labeling and R IUPAC ambiguous codes.
- `${outdir}/3_Aligned-Labelled-Aligned/` — Reads that have been labelled and undergone a second alignment step. Also includes original on/off target reads.
- `${outdir}/4_Final-bams/` — The final bams that have undergone all analysis steps.
- `${outdir}/5_DAF-QC/` — Transition stats for each sample.
- `${outdir}/6_pileups/` — Per sample bigwig files (only if `--create_bigwigs`).

---


<h3>Example runs</h3>

Local run (Docker), align + MSP/QC (minimum required parameters)

```bash
nextflow run main.nf \
  --sample_sheet inputs/sample_sheet.tsv \
  --ref_sheet inputs/reference_sheet.tsv \
  --target_sheet inputs/target_sheet.tsv \
  -profile local
```

Local run with bigwigs

```bash
nextflow run main.nf \
  --sample_sheet inputs/sample_sheet.tsv \
  --ref_sheet_path inputs/reference_sheet.tsv \
  --target_sheet inputs/target_sheet.tsv \
  --create_bigwigs true \
  --outdir results/full_run \
  -profile local
```

---

Tips

- Use absolute paths in `inputs/reference_sheet.tsv` 
- If you enable `-profile debug`, Nextflow will emit
  trace/timeline/report/dag files; use `--debug true` if you also want the
  channel `view()` messages.

