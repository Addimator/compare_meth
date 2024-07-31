library(GenomicRanges)
library(ChIPseeker)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

metilene <- read.table(snakemake@input[['metilene']], header=FALSE, sep="\t")
colnames(metilene) <- c("chr", "start", "end", "q_value", "mean_methylation_difference", "num_CpGs", "p_MWU", "p_2D_KS", "mean_g1", "mean_g2")
metilene <- metilene[order(metilene$chr, metilene$start), ]

# Rename chromosomes to be able to compare them to human TxDb
metilene$chr <- paste0("chr", metilene$chr)

# Create GRanges Object from DMRs
gr <- GRanges(seqnames = Rle(metilene$chr),
              ranges = IRanges(start = metilene$start, end = metilene$end),
              strand = Rle("*"))



# Annotate with chipseeker
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# print for debugging because not all the chromosomes are the same
unique(metilene$chr)
print("")
seqlevels(txdb)

chipseeker_output <- annotatePeak(gr, TxDb = txdb)
chipseeker_output_df <- as.data.frame(chipseeker_output)

# Create filtered output
chipseeker_output_df <- chipseeker_output_df %>%
  rename(chr = seqnames)
output <- merge(chipseeker_output_df, metilene, by = c("chr", "start", "end"))
output <- output %>%
  dplyr::select(chr, start, end, mean_methylation_difference, geneChr, geneStart, geneEnd)

# Speichern der annotierten Peaks als TSV-Datei
write.table(output, file = snakemake@output[['chipseeker']], sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)