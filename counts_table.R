library(Rsubread)
setwd("~/Desktop/BT549/8_mergedBAM")
samples = c("KH9_S1.sort.markdup.readgrps.merge.bam",
            "KH10_S2.sort.markdup.readgrps.merge.bam",
            "KH11_S3.sort.markdup.readgrps.merge.bam",
            "KH12_S4.sort.markdup.readgrps.merge.bam",
            "KH13_S5.sort.markdup.readgrps.merge.bam",
            "KH14_S6.sort.markdup.readgrps.merge.bam",
            "KH15_S7.sort.markdup.readgrps.merge.bam",
            "KH16_S8.sort.markdup.readgrps.merge.bam")
fc <- featureCounts(files=samples, 
                    annot.ext = 
                      "/home/mjaffe7/Desktop/BT549/3_genome/annotation.gtf",
                    isGTFAnnotationFile = TRUE,
                    isPairedEnd = FALSE,
                    genome = 
                      "/home/mjaffe7/Desktop/BT549/3_genome/GRCh38.primary_assembly.genome.fa")
rownames(fc$counts) <- gsub("\\..*", "", rownames(fc$counts))
colnames(fc$counts) <- gsub(".sort.markdup.readgrps.merge.bam", "", colnames(fc$counts))
save(fc, file = "fc.Rdata")
