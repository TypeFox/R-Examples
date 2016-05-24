

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Load prerequisites
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

require(refGenome)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Initialize example data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

ensfile <- system.file("extdata",
                       "hs.ensembl.62.small.RData",
                       package = "refGenome", mustWork=TRUE)

# A) Load Ensembl genome
ens <- loadGenome(ensfile)
enex<-refExons(ens)
gp <- getGenePositions(ens)
junc<-getSpliceTable(ens)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Run tests
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
source("test_overlap_juncs.r")
source("test_read_gtf.r")
source("test_unify_ranges.r")
source("test_unify_juncs.r")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Cleanup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

rm(ensfile, ens, gp, junc)
gc()

cat("[refGenome] rest-all.R tests finished.\n")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# END OF FILE
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
