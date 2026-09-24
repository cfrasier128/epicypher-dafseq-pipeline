process get_transition_stats {
    publishDir "$params.outdir/5_DAF-QC", mode: 'copy'
    cpus 1
    memory '4GB'
    container 'cfrasier/epi-dafseq:latest'

    input:
    tuple val(sample_id), path(bam_file), val(ref_name), val(target_name), val(target_precision), path(bam_index)
    output:
    tuple path("${sample_id}.transition_stats.${target_precision}.txt"), path("${sample_id}.read_transition.${target_precision}.tsv")
    script:
    """
    samtools view $bam_file | awk 'BEGIN{OFS="\\t"} {
      for(i=12;i<=NF;i++){
        if (\$i ~ /^CT:/) { split(\$i,a,":"); print \$1,a[3]; next }
        else if (\$i ~ /^GA:/) { split(\$i,a,":"); print \$1,a[3]; next }
      }
    }' > ${sample_id}.read_transition.${target_precision}.tsv

    awk '{sum+=\$2; n++} END { if (n>0) printf "%.6f\\n", sum/n; else print "NA" }' ${sample_id}.read_transition.${target_precision}.tsv > ${sample_id}.transition_stats.${target_precision}.txt
    """
}


process sort_index_bams {
    publishDir "$params.outdir/4_Final-bams/${target_precision}", mode: 'copy'
    cpus 4
    memory '8GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(sample_id), path(bam_file), val(ref_name), val(target_name), val(target_precision)
    output:
    tuple val(sample_id), path("${sample_id}.${target_precision}.second_alignment.aligned.nucs.sorted.bam"), val(ref_name), val(target_name), val(target_precision), path("${sample_id}.${target_precision}.second_alignment.aligned.nucs.sorted.bam.bai")
    script:
    """
    samtools view -F 2048 -b -o temp.bam $bam_file;
    samtools sort temp.bam -o ${sample_id}.${target_precision}.second_alignment.aligned.nucs.sorted.bam;
    samtools index -@ 8 ${sample_id}.${target_precision}.second_alignment.aligned.nucs.sorted.bam    
    """
}