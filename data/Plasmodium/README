# This directory contains RNA-Seq data (unpublished):
#   Alba3 KO dataset 
#   Done with Plasmodium berghei 
# Samples: 
#  820LT – wildtype 16hrs post invasion  
#  820Sch – wildtype 22-24 hrs post invasion 
#  G260 - Alba3 knock out 16 hours post invasion 
#  260Sch - Alba3 knock out 22-24 hours post invasion 
#  G260sch - Alba3 knock out 36 hours post invasion 
#  401Sch - non-gametocyte producer line which you can remove 

   mkdir -p data/Plasmodium/raw

# QC:
# These have been QC filtered with trim-galore

  mkdir -p analysis/Plasmodium/qc/
  for accession in $(cat data/Plasmodium/list_ids.txt)
  do
    trim_galore \
      data/Plasmodium/raw/$accession*.fastq.gz \
      --paired \
      --output_dir analysis/Plasmodium/qc/ \
      --basename $accession \
      --fastqc \
      --no_report_file
  done

# References:
# To perform RNA-Seq analysis we need reference sequences and annotations for each strain.
#
# For P. berghei, the reference can be found here:
# https://plasmodb.org/common/downloads/release-68/PbergheiANKA/
# and is associated with the following publication:
#	Fougère et al., Variant Exported Blood-Stage Proteins Encoded by Plasmodium Multigene Families Are Expressed in Liver Stages 
#            Where They Are Exported into the Parasitophorous Vacuole. PLoS Pathog 2016;12(11):e1005917
# We downloaded the following on 20/05/2024:
# https://plasmodb.org/common/downloads/release-68/PbergheiANKA/gff/data/PlasmoDB-68_PbergheiANKA.gff
# https://plasmodb.org/common/downloads/release-68/PbergheiANKA/fasta/data/PlasmoDB-68_PbergheiANKA_Genome.fasta
# https://plasmodb.org/common/downloads/release-68/PbergheiANKA/fasta/data/PlasmoDB-68_PbergheiANKA_AnnotatedTranscripts.fasta
#
# These have been stored in gzip format but can be unzipped using the following

   gunzip -c data/Plasmodium/reference/PlasmoDB-68_PbergheiANKA.gff.gz > data/Plasmodium/reference/PlasmoDB-68_PbergheiANKA.gff
   gunzip -c data/Plasmodium/reference/PlasmoDB-68_PbergheiANKA_Genome.fasta.gz > data/Plasmodium/reference/PlasmoDB-68_PbergheiANKA_Genome.fasta
   gunzip -c data/Plasmodium/reference/PlasmoDB-68_PbergheiANKA_AnnotatedTranscripts.fasta.gz > data/Plasmodium/reference/PlasmoDB-68_PbergheiANKA_AnnotatedTranscripts.fasta

# STAR alignment:
# As per the example notebook, we first use STAR to index the reference sequences and then run on the QC-filtered RNA-Seq FASTA files

  mkdir -p analysis/Plasmodium/star/ref
  STAR --runThreadN 4 \
    --runMode genomeGenerate \
    --genomeDir analysis/Plasmodium/star/ref \
    --genomeFastaFiles data/Plasmodium/reference/*_Genome.fasta \
    --sjdbGTFfile data/Plasmodium/reference/*.gff \
    --sjdbOverhang 75 \
    --genomeSAindexNbases 11

  for accession in $(cat data/Plasmodium/list_ids.txt)
  do
    mkdir -p analysis/Plasmodium/star/$accession
    STAR \
      --genomeDir analysis/Plasmodium/star/ref \
      --runThreadN 4 \
      --readFilesIn <(gunzip -c analysis/Plasmodium/qc/$accession*.fq.gz) \
      --outFileNamePrefix analysis/Plasmodium/star/$accession/$accession \
      --outSAMtype BAM SortedByCoordinate \
      --limitBAMsortRAM 1200000000 \
      --outSAMattributes Standard \
      --quantMode TranscriptomeSAM GeneCounts
  done
