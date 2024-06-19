# This directory contains data from the following publication:
# 	Protasio AV, Tsai IJ, Babbage A, Nichol S, Hunt M, Aslett MA, De Silva N, Velarde GS,
# Anderson TJC, Clark RC, et al. 2012. A systematically improved high quality genome
# and transcriptome of the human blood fluke Schistosoma mansoni. PLoS Negl Trop Dis 6: e1455.
#
# Full metadata for this experiment can be found at https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-451/sdrf and 
# https://www.ebi.ac.uk/ena/browser/view/PRJEB2350?show=related-records
#
# The full raw RNA-Seq data was downloaded from ENA Project PRJEB2350
# and downsampled to 1m reads in each file using rasusa in order to be
# tractable for the initial analysis within a limited resource environment.
# 
# Download:
  cd data/Schistosoma_mansoni
  mkdir -p raw
  cd raw

  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022880/ERR022880_1.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022874/ERR022874_2.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022877/ERR022877_1.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022873/ERR022873_1.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022878/ERR022878_2.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022881/ERR022881_1.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022873/ERR022873_2.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022878/ERR022878_1.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022872/ERR022872_1.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022875/ERR022875_1.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022879/ERR022879_2.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022876/ERR022876_2.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022882/ERR022882_2.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022883/ERR022883_1.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022879/ERR022879_1.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022877/ERR022877_2.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022880/ERR022880_2.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022872/ERR022872_2.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022882/ERR022882_1.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022876/ERR022876_1.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022875/ERR022875_2.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022883/ERR022883_2.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022881/ERR022881_2.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022874/ERR022874_1.fastq.gz
  cd ..

# Subsample
  mkdir -p subsampled
  for i in 72 73 74 75 76 77 78 79 80 81 82 83; do rasusa -n 1m -i raw/ERR0228"$i"_*.fastq.gz -o subsampled/ERR0228"$i"_1.fastq.gz subsampled/ERR0228"$i"_2.fastq.gz; done
  cd ../..


# The list of accessions for the files are saved to a file list_ids.txt

# QC:
# These have been QC filtered with trim-galore

  mkdir -p analysis/Schistosoma_mansoni/qc/
  for accession in $(cat data/Schistosoma_mansoni/list_ids.txt)
  do
    trim_galore \
      data/Schistosoma_mansoni/raw/$accession*.fastq.gz \
      --paired \
      --output_dir analysis/Schistosoma_mansoni/qc/ \
      --basename $accession \
      --fastqc \
      --no_report_file
  done

# References:
# To perform RNA-Seq analysis we need reference sequences and annotations for each strain.
#
# For S. mansoni, the reference can be found here:
# https://parasite.wormbase.org/Schistosoma_mansoni_prjea36577/
# and is associated with the following publications:
#       Berriman M et al., The genome of the blood fluke Schistosoma mansoni. Nature, 2009;460(7253):352-358
#       Protasio AV, et al., A systematically improved high quality genome and transcriptome of the human blood fluke Schistosoma mansoni. PLoS Negl Trop Dis, 2012;6(1):e1455
#            Where They Are Exported into the Parasitophorous Vacuole. PLoS Pathog 2016;12(11):e1005917
# We downloaded the following on 20/05/2024:
# https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS19.genomic.fa.gz
# https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3.gz
# https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS19.mRNA_transcripts.fa.gz
#
# These have been stored in gzip format but can be unzipped using the following

   gunzip -c data/Schistosoma_mansoni/reference/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3.gz > data/Schistosoma_mansoni/reference/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3
   gunzip -c data/Schistosoma_mansoni/reference/schistosoma_mansoni.PRJEA36577.WBPS19.genomic.fa.gz > data/Schistosoma_mansoni/reference/schistosoma_mansoni.PRJEA36577.WBPS19.genomic.fa
   gunzip -c data/Schistosoma_mansoni/reference/schistosoma_mansoni.PRJEA36577.WBPS19.mRNA_transcripts.fa.gz > data/Schistosoma_mansoni/reference/schistosoma_mansoni.PRJEA36577.WBPS19.mRNA_transcripts.fa

# The gff3 format file seems to have too nested complexity to provide gene level counts correctly so safer to convert to gtf first
   gffread -E data/Schistosoma_mansoni/reference/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3 -T -o- | more > data/Schistosoma_mansoni/reference/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gtf

# STAR alignment:
# As per the example notebook, we first use STAR to index the reference sequences and then run on the QC-filtered RNA-Seq FASTA files

  mkdir -p analysis/Schistosoma_mansoni/star/ref
  STAR --runThreadN 4 \
    --runMode genomeGenerate \
    --genomeDir analysis/Schistosoma_mansoni/star/ref \
    --genomeFastaFiles data/Schistosoma_mansoni/reference/*genomic.fa \
    --sjdbGTFfile data/Schistosoma_mansoni/reference/*.gtf \
    --sjdbOverhang 75 \
    --genomeSAindexNbases 11

  for accession in $(cat data/Schistosoma_mansoni/list_ids.txt)
  do
    mkdir -p analysis/Schistosoma_mansoni/star/$accession
    STAR \
      --genomeDir analysis/Schistosoma_mansoni/star/ref \
      --runThreadN 4 \
      --readFilesIn <(gunzip -c analysis/Schistosoma_mansoni/qc/$accession*.fq.gz) \
      --outFileNamePrefix analysis/Schistosoma_mansoni/star/$accession/$accession \
      --outSAMtype BAM SortedByCoordinate \
      --limitBAMsortRAM 1200000000 \
      --outSAMattributes Standard \
      --quantMode TranscriptomeSAM GeneCounts
  done

# MultiQC
# We combined the log files from STAR using multiqc with the following command

multiqc --outdir analysis/Schistosoma_mansoni/multiqc analysis/Schistosoma_mansoni/star/
