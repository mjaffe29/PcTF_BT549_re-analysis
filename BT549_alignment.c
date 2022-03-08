// To download off of GCP - 
//install GSUtil
  conda create --name
  conda activate name
  conda install -c conda-forge gsutil
  
// Configure GSUtil
gustil config

//Follow instructions - go to the link and activate account and then co to console and copy project id
//Then copy and paste the command that google cloud gives you into unix
  
// Sample names and corresponding samples:
// KH9_S1 —> untreated rep 1
// KH10_S2 —> untreated rep 2
// KH11_S3 —> treated 24hrs rep 1
// KH12_S4 —> treated 24hrs rep 2
// KH13_S5 —> treated 48hrs rep 1
// KH14_S6 —> treated 48hrs rep 2
// KH15_S7 —> treated 72hrs rep 1
// KH16_S8 —> treated 72hrs rep 2

//trimming
  parallel -j20 java -jar /home/mjaffe7/Desktop/RNA-seq_tools/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 {1}_L001_R1_001.fastq /home/mjaffe7/Desktop/BT549/trimmed/{1}_L001_Trim_430minlen50.fastq ILLUMINACLIP:/home/mjaffe7/Desktop/RNA-seq_tools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:25 SLIDINGWINDOW:4:30 MINLEN:50 ::: KH9_S1 KH10_S2 KH11_S3 KH12_S4 KH13_S5 KH14_S6 KH15_S7 KH16_S8; parallel -j20 java -jar /home/mjaffe7/Desktop/RNA-seq_tools/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 {1}_L002_R1_001.fastq /home/mjaffe7/Desktop/BT549/trimmed/{1}_L002_Trim_430minlen50.fastq ILLUMINACLIP:/home/mjaffe7/Desktop/RNA-seq_tools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:25 SLIDINGWINDOW:4:30 MINLEN:50 ::: KH9_S1 KH10_S2 KH11_S3 KH12_S4 KH13_S5 KH14_S6 KH15_S7 KH16_S8; parallel -j20 java -jar /home/mjaffe7/Desktop/RNA-seq_tools/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 {1}_L003_R1_001.fastq /home/mjaffe7/Desktop/BT549/trimmed/{1}_L003_Trim_430minlen50.fastq ILLUMINACLIP:/home/mjaffe7/Desktop/RNA-seq_tools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:25 SLIDINGWINDOW:4:30 MINLEN:50 ::: KH9_S1 KH10_S2 KH11_S3 KH12_S4 KH13_S5 KH14_S6 KH15_S7 KH16_S8; parallel -j20 java -jar /home/mjaffe7/Desktop/RNA-seq_tools/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 {1}_L004_R1_001.fastq /home/mjaffe7/Desktop/BT549/trimmed/{1}_L004_Trim_430minlen50.fastq ILLUMINACLIP:/home/mjaffe7/Desktop/RNA-seq_tools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:25 SLIDINGWINDOW:4:30 MINLEN:50 ::: KH9_S1 KH10_S2 KH11_S3 KH12_S4 KH13_S5 KH14_S6 KH15_S7 KH16_S8

//fastqc
  conda install -c bioconda fastqc
  for FILE in *; do fastqc $FILE; done
//multiqc      
  conda install -c bioconda multiqc
//will need this following line:
  export PYTHONIOENCODING=utf-8
  multiqc .
    
//install necessary tools
    conda install -c bioconda star
    conda install -c bioconda samtools
     
//download most recent genome from GENCODE using wget and the ftp link
//running STAR
    STAR --genomeDir /home/mjaffe7/Desktop/BT549/3_genome --runMode genomeGenerate --genomeFastaFiles /home/mjaffe7/Desktop/BT549/3_genome/GRCh38.primary_assembly.genome.fa --sjdbGTFfile annotation.gtf --runThreadN 12
    
    picard CreateSequenceDictionary R=GRCh38.primary_assembly.genome.fa O=GRCh38.primary_assembly.genome.fa.dict

//ALL OF THE FILES THE WHILE LOOPS ARE BASED OFF OF ARE IN THE GITHUB
//aligning —> taking the aligned reads
    while IFS= read -r line; do 
      STAR --genomeDir /home/mjaffe7/Desktop/BT549/3_genome --sjdbGTFfile /home/mjaffe7/Desktop/BT549/3_genome/annotation.gtf --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --readFilesIn /home/mjaffe7/Desktop/BT549/2_trimmed/${line}.fastq --outSAMtype BAM Unsorted --outFileNamePrefix /home/mjaffe7/Desktop/BT549/4_STAR_out/${line}_pass1. —runThreadN 8; 
    done < files.txt
          
//sorting bam files
    while IFS= read -r line; do bamtools sort -in /home/mjaffe7/Desktop/BT549/4_STAR_out/${line}_Trim_430minlen50_pass1.Aligned.out.bam -out /home/mjaffe7/Desktop/BT549/5_BAM/${line}.bam
        picard MarkDuplicates INPUT=/home/mjaffe7/Desktop/BT549/5_BAM/${line}.bam OUTPUT=/home/mjaffe7/Desktop/BT549/5_BAM/${line}.markdup.bam METRICS_FILE=/home/mjaffe7/Desktop/BT549/5_BAM/${line}.markdup.picardMetrics.txt REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT; 
    done < files_bam.txt

//marking duplicates
    while IFS= read -r line; 
      do picard MarkDuplicates INPUT=/home/mjaffe7/Desktop/BT549/5_BAM/${line}.bam OUTPUT=${line}.sort.markdup.bam METRICS_FILE=/home/mjaffe7/Desktop/BT549/5_BAM/${line}.markdup.picardMetrics.txt REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT; 
    done < files_bam.txt

//adding read groups
    while IFS=, read -r key value; 
      do picard AddOrReplaceReadGroups INPUT=/home/mjaffe7/Desktop/BT549/6_sorted_BAM/${key}.sort.markdup.bam OUTPUT=${key}.sort.markdup.readgrps.bam RGID=${key} RGSM=${value} RGPL=ILLUMINA RGPU=lib1 RGLB=${key} RGCN=ArizonaStateUniversity RGDS=homosapien VALIDATION_STRINGENCY=LENIENT; 
    done < cell_type.txt

//obtaining stats, see BAMtools stats pt1 for the output  
    while IFS= read -r line; 
      do bamtools stats -in ${line}.sort.markdup.readgrps.bam; 
    done < files_bam.txt
                    
//merging the lanes
    while IFS=, read -r f1 f2 f3 f4 f5; 
      do bamtools merge -in /home/mjaffe7/Desktop/BT549/7_readgroups/${f1}.sort.markdup.readgrps.bam -in /home/mjaffe7/Desktop/BT549/7_readgroups/${f2}.sort.markdup.readgrps.bam -in /home/mjaffe7/Desktop/BT549/7_readgroups/${f3}.sort.markdup.readgrps.bam -in /home/mjaffe7/Desktop/BT549/7_readgroups/${f4}.sort.markdup.readgrps.bam -out ${f5}.sort.markdup.readgrps.merge.bam; 
    done < files_bam_2.txt

//creatin indexes such that they can be read by IGV                                                                                                                                                        
    while IFS= read -r line; 
      do bamtools index -in ${line}.sort.markdup.readgrps.merge.bam; 
    done < files_bam_3.txt

//obtaining stats, see BAMtools stats pt2 for output                                                                                                                                                        
    while IFS= read -r line; 
      do bamtools stats -in ${line}.sort.markdup.readgrps.merge.bam; 
    done < files_bam_3.txt
