nextflow.enable.dsl = 2


process DEEP_VARIANT {
    container "google/deepvariant:1.0.0"
    cpus 16
    memory "32 GB"


    input:
    path(ref_fasta_gz)
    path(ref_fasta_fai)
    path(chr_bam)
    path(chr_bam_bai)

    output:
    path("output")

    script:

   """
   mkdir ./output
   
   /opt/deepvariant/bin/run_deepvariant \
      --model_type WGS \
      --ref GRCh38_no_alt_analysis_set.fasta \
      --reads HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam \
      --output_vcf ./output/HG002.output.vcf.gz \
      --output_gvcf ./output/HG002.output.g.vcf.gz \
      --num_shards ${task.cpus} \
      --regions chr20 
      
   """

}


workflow wgs_test {
// References
    ref_fasta = "s3://nf-core-awsmegatests/deepvariant/test_data/wgs/reference/GRCh38_no_alt_analysis_set.fasta"
    ref_fasta_fai = "s3://nf-core-awsmegatests/deepvariant/test_data/wgs/reference/GRCh38_no_alt_analysis_set.fasta.fai"

// HG002 chr20 BAM
    chr_bam = "s3://nf-core-awsmegatests/deepvariant/test_data/wgs/input/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam"
    chr_bam_bai = "s3://nf-core-awsmegatests/deepvariant/test_data/wgs/input/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam.bai"


    DEEP_VARIANT(ref_fasta, ref_fasta_fai, chr_bam, chr_bam_bai)
}
