library(trackObjs)
source("trackObjs/inst/performanceTrials/funsForTesting.R")
options(width=130)
sessionInfo()
gc(reset=TRUE)
# Create a tracking directory with 33 large objects in
runSaveLoadTest("create", scale=12, nObjs=33, clobber=TRUE, options=list(cache=FALSE))
gc()
gc(reset=TRUE)
# Retrieve all objects from the tracking directory, but don't keep objects in memory
runSaveLoadTest("verify", options=list(cache=FALSE))
gc()
gc(reset=TRUE)
# Retrieve all objects from the tracking directory, and do keep objects in memory
# Note that after this runs, memory consumption is much higher.
runSaveLoadTest("verify", options=list(cache=TRUE))
gc()
# Remove the tracking dir
unlink("test1", recursive=TRUE)
