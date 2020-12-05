nextflow.enable.dsl = 2

// References
REF_FTPDIR = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids"
ref_fasta_gz = "${REF_FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
ref_fasta_fai = "${REF_FTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai"

// Genomes
READS_FTPDIR = "ftp://ftp.ncbi.nlm.nih.gov//giab/ftp/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020"
bed_file = "${READS_FTPDIR}/HG002_GRCh38_1_22_v4.2_benchmark.bed"
vcf_file = "${READS_FTPDIR}/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz"
vcf_index_file = "${READS_FTPDIR}/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi"

// HG002 chr20 BAM
HTTPDIR = "https://storage.googleapis.com/deepvariant/case-study-testdata"
chr_bam = "${HTTPDIR}/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam"
chr_bam_bai = "${HTTPDIR}/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam.bai"



process DOWNLOAD_FILES {
    input:
    file(ref_fasta_gz)
    file(ref_fasta_fai)
    file(bed_file)
    file(vcf_file)
    file(vcf_index_file)
    file(chr_bam)
    file(chr_bam_bai)

    output:
    path("reference")
    path("benchmark")
    path("input")



    script:
    """

    mkdir reference
    mv $ref_fasta_gz reference
    mv $ref_fasta_fai reference

    mkdir benchmark
    mv $bed_file benchmark
    mv $vcf_file benchmark
    mv $vcf_index_file  benchmark
    
    mkdir input
    mv $chr_bam input
    mv $chr_bam_bai input

    """

}


process DEEP_VARIANT {
    container "google/deepvariant:1.0.0"

    input:
    path(reference)
    path(benchmark)
    path(input)

    output:
    path("*vcf.gz")

    shell:

    '''
    mkdir output
   
   /opt/deepvariant/bin/run_deepvariant \
      --model_type WGS \
      --ref /reference/GRCh38_no_alt_analysis_set.fasta \
      --reads /input/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam \
      --output_vcf /output/HG002.output.vcf.gz \
      --output_gvcf /output/HG002.output.g.vcf.gz \
      --num_shards $(nproc) \
      --regions chr20 
      
   '''

}


workflow {
    DOWNLOAD_FILES()
    DEEP_VARIANT(DOWNLOAD_FILES.out)
}
