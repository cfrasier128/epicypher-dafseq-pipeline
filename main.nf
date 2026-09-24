#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// inputs
params.sample_sheet = ''
params.ref_sheet= ''
params.outdir = "${workflow.launchDir}/results"
params.confidence_ml_val = '250'
params.minimum_msp_dist = '10'
params.target_sheet = ''

params.split_reads = false
params.create_bigwigs = false

include { create_bigwigs } from './subworkflows/create_bigwigs.nf'
include { get_transition_stats as get_transition_stats_on_target } from './modules/epicypher-dafseq-processes.nf'
include { get_transition_stats as get_transition_stats_off_target } from './modules/epicypher-dafseq-processes.nf'
include { sort_index_bams as sort_index_bams_on_target } from './modules/epicypher-dafseq-processes.nf'
include { sort_index_bams as sort_index_bams_off_target } from './modules/epicypher-dafseq-processes.nf'

process align_reads{
    publishDir "$params.outdir/1_Aligned-NotLabelled", mode: 'copy'
    cpus 8
    memory '16GB'
    container 'cfrasier/epi-dafseq:latest'

    input:
    tuple val(sample_id), path(fastq), val(ref_name), val(target_name), path(ref_fasta)
    output:
    tuple val(sample_id), path("${sample_id}.daf.aligned.bam"), val(ref_name), val(target_name)
    script:
    """
    minimap2 -d ${ref_name}.mmi $ref_fasta;
    minimap2 --MD -Y -y -a -x map-ont ${ref_name}.mmi $fastq | samtools view -b > ${sample_id}.daf.aligned.bam
    """
}

process label_reads{
    publishDir "$params.outdir/2_Aligned-Labelled", mode: 'copy'
    cpus 8
    memory '16GB'
    container 'cfrasier/epi-dafseq:latest'

    input:
    tuple val(sample_id), path(bam_file), val(ref_name), val(target_name)
    output:
    tuple val(sample_id), path("${sample_id}.daf.labelled.bam"), val(ref_name), val(target_name)
    script:
    """
    conda run -n dafseq python3 /opt/dafseq/labelReads.py -i $bam_file -o ${sample_id}.daf.labelled.bam
    """
}   

process align_labelled_reads{
    publishDir "$params.outdir/3_Aligned-Labelled-Aligned", mode: 'copy'
    cpus 8
    memory '16GB'
    container 'cfrasier/epi-dafseq:latest'

    input:
    tuple val(sample_id), path(labelled_bam), val(ref_name), val(target_name) , path(ref_fasta), path(ref_fai)
    output:
    tuple val(sample_id), path("${sample_id}.no_MM_ML.bam"), val(ref_name), val(target_name)
    script:
    """
    minimap2 -d ${ref_name}.mmi $ref_fasta;
    minimap2 --MD -Y -a -y -x map-pb ${ref_name}.mmi <(samtools fastq -T '*' $labelled_bam) | samtools view -bh -x "MM,ML" - > ${sample_id}.no_MM_ML.bam    
    """
}

process convert_6ma {
    publishDir "$params.outdir/3_Aligned-Labelled-Aligned", mode: 'copy'
    cpus 8
    memory '16GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(sample_id), path(bam_file), val(ref_name), val(target_name)
    output:
    tuple val(sample_id), path("${sample_id}.no_MM_ML.ddda.bam"), val(ref_name), val(target_name)
    script:
    """
    ft ddda-to-m6a $bam_file ${sample_id}.no_MM_ML.ddda.bam
    """
}   

process get_on_target_only {
    publishDir "$params.outdir/3_Aligned-Labelled-Aligned", mode: 'copy'
    cpus 8
    memory '16GB'
    container 'cfrasier/epi-dafseq:0.20'

    input:
    tuple val(sample_id), path(bam_file), val(ref_name), val(target_name), val(target_chr), val(target_start), val(target_end)
    output:
    tuple val(sample_id), path("${sample_id}.on_target_only.bam"), val(ref_name), val(target_name), val("on_target"), emit: on_target_bam
    tuple val(sample_id), path("${sample_id}.off_target_only.bam"), val(ref_name), val(target_name), val("off_target"), emit: off_target_bam
    script:
    """
    echo -e "${target_chr}\t${target_start}\t${target_end}" > ${sample_id}.on_target.bed
    bedtools intersect -abam $bam_file -b ${sample_id}.on_target.bed > ${sample_id}.on_target_only.bam
    bedtools intersect -abam $bam_file -b ${sample_id}.on_target.bed -v > ${sample_id}.off_target_only.bam
    """
}

process add_nucleosomes {
    publishDir "$params.outdir/3_Aligned-Labelled-Aligned", mode: 'copy'
    cpus 8
    memory '16GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(sample_id), path(bam_file), val(ref_name), val(target_name)
    output:
    tuple val(sample_id), path("${sample_id}.second_alignment.aligned.nucs.bam"), val(ref_name), val(target_name)
    script:
    """
    ft add-nucleosomes -n 60 -c 70 --min-distance-added 15 -d 10 $bam_file > ${sample_id}.second_alignment.aligned.nucs.bam;
    """
} 


process split_reads {
    publishDir "$params.outdir/4_Final-bams", mode: 'copy'
    cpus 4
    memory '8GB'
    container 'cfrasier/epi-dafseq:latest'

    input:
    tuple val(sample_id), path(bam_file), val(ref_name), val(target_name), val(target_precision), path(bam_index)
    output:
    tuple val("${sample_id}_CT"), path("${sample_id}_CT.second_alignment.aligned.nucs.sorted.split.bam"), val(ref_name), path("${sample_id}_CT.second_alignment.aligned.nucs.sorted.split.bam.bai"), emit: ct_split_bam
    tuple val("${sample_id}_GA"), path("${sample_id}_GA.second_alignment.aligned.nucs.sorted.split.bam"), val(ref_name), path("${sample_id}_GA.second_alignment.aligned.nucs.sorted.split.bam.bai"), emit: ga_split_bam
    script:
    """
    samtools view -h $bam_file | grep -Pe "\tCT:|^@" | samtools view -b -o ${sample_id}_CT.second_alignment.aligned.nucs.sorted.split.bam -;
    samtools view -h $bam_file | grep -Pe "\tGA:|^@" | samtools view -b -o ${sample_id}_GA.second_alignment.aligned.nucs.sorted.split.bam -;
    samtools index -@ 8 ${sample_id}_CT.second_alignment.aligned.nucs.sorted.split.bam;
    samtools index -@ 8 ${sample_id}_GA.second_alignment.aligned.nucs.sorted.split.bam
    """
}



workflow{

    target_sheet_ch = channel.fromPath("${params.target_sheet}")
        .splitCsv(skip: 1, sep: '\t')
        .map { row -> tuple(row[0], row[1], row[2], row[3]) }

    // references_ch layout (post map) -> ref_fasta, ref_fai, ref_name
    references_ch = channel.fromPath("${params.ref_sheet}")
        .splitCsv(skip: 1, sep: '\t')
        .map { row -> tuple(row[1], row[2], row[0]) }

    sample_sheet_ch = channel.fromPath("${params.sample_sheet}")
        .splitCsv(skip: 1, sep: '\t')
        .map { row -> tuple(row[0], file(row[1]), row[2], row[3]) }
    // sample_sheet_ch layout (post map) -> samp_name, fastq_paths, ref_name, target_name

    aligned_bams_input_ch = sample_sheet_ch
        .combine(references_ch, by: 2)
        .map { row -> tuple("${row[1]}.${row[3]}", row[2], row[0], row[3], row[4]) }
    // -> samp_name, fastq_paths, ref_name, target_name, ref_fasta

    align_reads(aligned_bams_input_ch)
    label_reads(align_reads.out)
    label_reads.out.combine(references_ch, by: 2).map { row -> tuple(row[1], row[2], row[0], row[3], row[4], row[5]) }.set { align_labelled_reads_input_ch }
    align_labelled_reads(align_labelled_reads_input_ch)
    convert_6ma(align_labelled_reads.out)



    add_nucleosomes(convert_6ma.out)
    add_nucleosomes.out
        .combine(target_sheet_ch.map { row -> tuple(row[0], row[1], row[2], row[3]) }, by: 3)
        .map { row -> tuple(row[1], row[2], row[3], row[0], row[4], row[5], row[6]) }
        .set { get_on_target_only_input_ch }
    // sampname, bam, ref_name, target_name, target_chr, target_start, target_end

    get_on_target_only(get_on_target_only_input_ch)
    sort_index_bams_on_target(get_on_target_only.out.on_target_bam)
    get_transition_stats_on_target(sort_index_bams_on_target.out)
    
    sort_index_bams_off_target(get_on_target_only.out.off_target_bam)
    get_transition_stats_off_target(sort_index_bams_off_target.out)

    if (params.split_reads) {
        split_reads(sort_index_bams_on_target.out)
        split_reads.out.ct_split_bam
                .mix(split_reads.out.ga_split_bam).set{ create_bigwigs_input_ch }
    }
    else {
        sort_index_bams_on_target.out.map{row -> tuple(row[0], row[1], row[2], row[5])}.set{ create_bigwigs_input_ch }
    }
    if (params.create_bigwigs) {
        create_bigwigs(create_bigwigs_input_ch.map { row -> tuple(row[0], row[1], row[2], row[3]) }, references_ch)
    }
}