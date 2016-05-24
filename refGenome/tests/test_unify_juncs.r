
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Requires: Initialized objects (as done by test-all.R header)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Artificially added transcript_biotype data from ensembl version 80
# For transcript_id's which are not found idn ensembl 80,
# the transcript_biotype is manually set to "nonsense_mediated_decay".
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# ensfile <- system.file("extdata", "hs.ensembl.62.small.RData",
#                                         package = "refGenome", mustWork=TRUE)
# ens <- loadGenome(ensfile)
# jtn <- ens@ev$gtf$transcript_id
# mtc <- match(jtn, en80@ev$gtf$transcript_id)
# biot <- en80@ev$gtf$transcript_biotype[mtc]
# biot[is.na(biot)] <- "nonsense_mediated_decay"
# ens@ev$gtf$transcript_biotype <- biot
# saveGenome(ens, 
#     file="~/projects/R/refGenome/inst/extdata/hs.ensembl.62.small.RData",
#     useBasedir=FALSE)
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


uj <- unifyJuncs(junc)

if(! all(uj@ev$gtf$cnNmd >= 0) )
    stop("[test_unify_juncs] All cnNmd must be >= 0")

if(! all(uj@ev$gtf$cnNmd <= uj@ev$gtf$nSites) )
    stop("[test_unify_juncs] All cnNmd must be <= nSites")

rm(uj)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# END OF FILE
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
