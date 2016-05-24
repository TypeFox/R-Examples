

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Load prerequisites
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

require(rbamtools)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
## Initialize example data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

bam<- system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx<- system.file("extdata", "accepted_hits.bam.bai", package="rbamtools")
# Open BAM-file for reading
reader<-bamReader(bam,idx=TRUE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Run tests
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source("test_bam_range.r")
source("test_bam_header.r")
source("test_bam_align.r")
source("test_range_seg_count.r")
source("test_genome_partition.r")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Cleanup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

bamClose(reader)
rm(reader)
gc()

cat("[rbamtools] rest-all.R tests finished.\n")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# END OF FILE
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
