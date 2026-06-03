#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_path = ''
params.input_string_filter = ''
params.outdir = params.input_path + '/dafseq_output'
params.reference_genome = 'T2T'
params.confidence_ml_val = '250'
params.minimum_msp_dist = '10'

process align_reads{
    publishDir "$params.outdir/1_Aligned-NotLabelled", mode: 'copy'
    cpus 8
    memory '16GB'
    container 'cfrasier/epi-dafseq:latest'

    input:
    tuple val(sample_id), path(fastq)
    path ref_mmi
    output:
    tuple val(sample_id), path("${sample_id}.daf.aligned.bam")
    script:
    """
    minimap2 --MD -Y -a -y -x map-ont $ref_mmi $fastq | samtools view -b > ${sample_id}.daf.aligned.bam
    """
}

process label_reads{
    publishDir "$params.outdir/2_Aligned-Labelled", mode: 'copy'
    cpus 8
    memory '16GB'
    container 'cfrasier/epi-dafseq:latest'

    input:
    tuple val(sample_id), path(bam_file)
    output:
    tuple val(sample_id), path("${sample_id}.daf.labelled.bam")
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
    tuple val(sample_id), path(labelled_bam)
    path ref_mmi
    output:
    tuple val(sample_id), path("${sample_id}.no_MM_ML.sam")
    script:
    """
    minimap2 --MD -Y -a -y -x map-pb $ref_mmi <(samtools fastq -T '*' $labelled_bam) | samtools view -h -x "MM,ML" - > ${sample_id}.no_MM_ML.sam    
    """
}

process convert_6ma {
    publishDir "$params.outdir/3_Aligned-Labelled-Aligned", mode: 'copy'
    cpus 8
    memory '16GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(sample_id), path(bam_file)
    output:
    tuple val(sample_id), path("${sample_id}.no_MM_ML.ddda.sam")
    script:
    """
    timeout --preserve-status 10 ft ddda-to-m6a $bam_file ${sample_id}.no_MM_ML.ddda.sam || true
    """
}   

process add_nucleosomes {
    publishDir "$params.outdir/3_Aligned-Labelled-Aligned", mode: 'copy'
    cpus 8
    memory '16GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(sample_id), path(bam_file)
    output:
    tuple val(sample_id), path("${sample_id}.ambiguous.aligned.nucs.bam")
    script:
    """
    ft add-nucleosomes -n 60 -c 70 --min-distance-added 15 -d 10 $bam_file > ${sample_id}.ambiguous.aligned.nucs.bam;
    """
} 

process sort_index_bams {
    publishDir "$params.outdir/4_Final-bams", mode: 'copy'
    cpus 4
    memory '8GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(sample_id), path(bam_file)
    output:
    tuple val(sample_id), path("${sample_id}.ambiguous.aligned.nucs.sorted.bam"), path("${sample_id}.ambiguous.aligned.nucs.sorted.bam.bai")
    script:
    """
    samtools sort $bam_file -o ${sample_id}.ambiguous.aligned.nucs.sorted.bam;
    samtools index -@ 8 ${sample_id}.ambiguous.aligned.nucs.sorted.bam    
    """
}

process get_transition_stats {
    publishDir "$params.outdir/5_DAF-QC", mode: 'copy'
    cpus 1
    memory '4GB'
    container 'cfrasier/epi-dafseq:latest'

    input:
    tuple val(sample_id), path(bam_file)
    output:
    tuple path("${sample_id}_transition_stats.txt"), path("${sample_id}_read_transition.tsv")
    script:
    """
    samtools view $bam_file | awk 'BEGIN{OFS="\\t"} {
      for(i=12;i<=NF;i++){
        if (\$i ~ /^CT:/) { split(\$i,a,":"); print \$1,a[3]; next }
        else if (\$i ~ /^GA:/) { split(\$i,a,":"); print \$1,a[3]; next }
      }
    }' > ${sample_id}_read_transition.tsv

    awk '{sum+=\$2; n++} END { if (n>0) printf "%.6f\\n", sum/n; else print "NA" }' ${sample_id}_read_transition.tsv > ${sample_id}_transition_stats.txt
    """
}

process create_pileups {
    publishDir "$params.outdir/4_Pileups_Bigwigs/1_Pileups/"
    cpus 1
    memory '4GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(sampname), path(aligned_bam), path(bam_index)
    output:
    tuple val(sampname), path("*.tsv"), emit: pileups
    script:
    """
    ft pileup \
    --m6a --ml $params.confidence_ml_val \
    --cpg \
    --per-base \
    -t $task.cpus \
    --ftx "len(msp)>$params.minimum_msp_dist" \
    $aligned_bam | awk -v OFS="\t" -v FS="\t" '{print \$1,\$2,\$3,\$9/(\$4+0.1),\$10/(\$4+0.1),\$7/(\$4+0.1)}' > ${sampname}.pileup_all.tsv
    """
}

process pileupbedgraphtobigwig_6ma{
    publishDir "$params.outdir/4_Pileups_Bigwigs/2_BigWigs/", mode: 'copy'
    cpus 1
    memory '4GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(sampname), path(pileup)
    path(genome_chromsizes)
    output:
    tuple val(sampname), path("*.bw")
    script:
    """
    cut -f 1,2,3,4 $pileup > temp.bedgraph
    bedGraphToBigWig temp.bedgraph $genome_chromsizes ${sampname}.perc6ma.bw
    """
}

process pileupbedgraphtobigwig_5mC{
    publishDir "$params.outdir/4_Pileups_Bigwigs/2_BigWigs/", mode: 'copy'
    cpus 1
    memory '4GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(sampname), path(pileup)
    path(genome_chromsizes)
    output:
    tuple val(sampname), path("*.bw")
    script:
    """
    cut -f 1-3,5 $pileup > temp.bedgraph
    bedGraphToBigWig temp.bedgraph $genome_chromsizes ${sampname}.perc5mc.bw
    """
}

process pileupbedgraphtobigwig_nuc{
    publishDir "$params.outdir/4_Pileups_Bigwigs/2_BigWigs/", mode: 'copy'
    cpus 1
    memory '4GB'
    container 'cfrasier/epi-fiberseq:latest'

    input:
    tuple val(sampname), path(pileup)
    path(genome_chromsizes)
    output:
    tuple val(sampname), path("*.bw")
    script:
    """
    cut -f 1-3,6 $pileup > temp.bedgraph
    bedGraphToBigWig temp.bedgraph $genome_chromsizes ${sampname}.percnuc.bw
    """
}

workflow{
    if (params.reference_genome == 'T2T') {
        ref_mmi = "/media/genomics/18Tb_1/references/T2T/T2T.ch13v2.0.mmi"
        ref_chromsizes = "/media/genomics/18Tb_1/references/T2T/hs1.chrom.sizes"
        ref_fasta = "/media/genomics/18Tb_1/references/T2T/chm13v2.0.clean.fasta"
        ref_fai = "/media/genomics/18Tb_1/references/T2T/chm13v2.0.clean.fasta.fai"
    }
    else if (params.reference_genome == 'hg38') {
        ref_mmi = "/media/genomics/18Tb_1/references/hg38/GCF_000001405.40_GRCh38.p14_genomic.mmi"
        ref_chromsizes = "/media/genomics/18Tb_1/references/hg38/hg38.uscsnames.v40.chrom.sizes"
        ref_fasta = "/media/genomics/18Tb_1/references/hg38/GCF_000001405.40_GRCh38.p14_genomic.fa"
        ref_fai = "/media/genomics/18Tb_1/references/hg38/GCF_000001405.40_GRCh38.p14_genomic.fa.fai"
    }
    else if (params.reference_genome == 'mm10') {
        ref_mmi = "/media/genomics/18Tb_1/references/mm10/GCF_000001635.27_GRCm39_genomic.mmi"
        ref_chromsizes = "/media/genomics/18Tb_1/references/mm10/GCF_000001635.27_GRCm39_genomic.chrom.sizes"
        ref_fasta = "/media/genomics/18Tb_1/references/mm10/GCF_000001635.27_GRCm39_genomic.fa"
        ref_fai = "/media/genomics/18Tb_1/references/mm10/GCF_000001635.27_GRCm39_genomic.fa.fai"
    }
    else if (params.reference_genome == 'LS_ref') {
        ref_mmi = "/media/genomics/18Tb_1/references/LS_dNuc/LS_ref.mmi"
        ref_chromsizes = "/media/genomics/18Tb_1/references/LS_dNuc/LS_ref.fasta.chrom.sizes"
        ref_fasta = "/media/genomics/18Tb_1/references/LS_dNuc/LS_ref.fasta"
        ref_fai = "/media/genomics/18Tb_1/references/LS_dNuc/LS_ref.fasta.fai"
    }
    else {
        exit 1
    }

    if (!params.input_path) {
        println "Please provide an input path with --input_path"
        exit 1
    }
    if (!params.input_string_filter) {
        input_ch = channel.fromPath("${params.input_path}/*.fastq")
    }
    else {
        input_ch = channel.fromPath("${params.input_path}/*${params.input_string_filter}*.fastq")
    }
    input_ch.map { file -> tuple( file.baseName.split('.fastq')[0], file ) }.set{input_fastq_names_ch}
    input_fastq_names_ch.view()
    align_reads(input_fastq_names_ch, ref_mmi)
    label_reads(align_reads.out)
    align_labelled_reads(label_reads.out, ref_mmi)
    convert_6ma(align_labelled_reads.out)
    add_nucleosomes(convert_6ma.out)
    sort_index_bams(add_nucleosomes.out)
    get_transition_stats(label_reads.out)
    create_pileups(sort_index_bams.out)
    pileupbedgraphtobigwig_6ma(create_pileups.out.pileups, ref_chromsizes)
    pileupbedgraphtobigwig_5mC(create_pileups.out.pileups, ref_chromsizes)
    pileupbedgraphtobigwig_nuc(create_pileups.out.pileups, ref_chromsizes)
}