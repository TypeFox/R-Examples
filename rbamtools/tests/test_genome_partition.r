
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Create genomePartition object
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Open (indexed) BAM file (done by calling test-all.R)
# bam<-system.file("extdata", "accepted_hits.bam", package="rbamtools")
# reader<-bamReader(bam,idx=TRUE)

# Provide exon positions
id <- 1:13
seqid <- "chr1"
gene <- "WASH7P"
ensg_id <- "ENSG00000227232"
start <- c(14411, 15000, 15796, 15904, 16607, 16748, 16858, 17233,
           17602, 17915, 18268, 24737, 29534)
end <-   c(14502, 15038, 15901, 15947, 16745, 16765, 17055, 17364,
           17742, 18061, 18366, 24891, 29806)

ref <- data.frame(id=id, seqid=seqid, begin=start, end=end, gene=gene, ensg=ensg_id)

# Create partition (adds equidistant grid)
partition <- genomePartition(reader, ref)

# data_frame test:
rn <- as.numeric(rownames(partition@ev$reflist[[1]]))
idx <- 1:length(rn)
if(any(rn!=idx))
    stop("[test_genome_partition.r] Error in rownames produced by data_frame!")

