# count-MAGE-seq

This readme file describes the  sequence data, bioinformatics steps and scripts used to analyze the MAGE-seq data in:

Konate, M., Plata, G., Park,J., Wang, H.H., Vitkup, D. Molecular function limits divergent protein evolution on planetary timescales.

Below we describe the steps used to obtain the number of reads associated with each mutant and the "wild type" strain for 
a particular sample. Each sample corresponds to either the first or the second half of the FolA protein, and a particular 
timepoint in the competition experiment.

Protocol:

1. Raw sequence reads can be accsessed from the SRA database with accession number: SRP152339 
   https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP152339
   
2. For each sample, paired reads were merged using SeqPrep-1.1 (https://github.com/jstjohn/SeqPrep) using the following command:
  
  SeqPrep -f ../T1_S1_L001_R1_001.fastq.gz -r ../T1_S1_L001_R2_001.fastq.gz -1 /dev/null -2 /dev/null \
  -s T1_S1_L001.merged.10.fastq.gz -L 100 -m 0.1
  
3. Merged reads were aligned to the corresponding reference sequence (first half or second half of FolA) using bowtie2
   (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml):
   
   bowtie2-build Reference_First_Half.fas Reference_First_Half
   bowtie2 -x Reference_First_Half -U T1_S1_L001.merged.10.fastq.gz -S T1.Merged.bowtie.10.local.sam -p 4 --local
   
4. Sequence positions that were identical to the reference were converted to = using the calmd command in samtools
   (http://samtools.sourceforge.net/):
   
   samtools calmd -e -S T1.Merged.bowtie.10.local.sam Reference_First_Half.fas > T1.Merged.bowtie.10.local.=.sam

5. A perl script was written to count the number of reads aligned to each single codon mutant and the "wild type":

  perl CountVariantsFirstHalf.pl T1.Merged.bowtie.10.local.=.sam T1.Merged.bowtie.10.local.sam Reference_First_Half.fas \
  T1.Merged.COUNTS.10.local.QUAL_13.txt QualitySymbols.txt 13
  

6. A similar procedure was used to analyze samples from the remaining timepoints and the second half of the protein, 
   changing the reference sequence and perl script as appropriate
   
