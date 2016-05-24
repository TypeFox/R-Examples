## --------------------------------------------------------------------------
## Test loading of saved rstream object
##

library(rstream)
load("test_load.RData")
rstream.packed(stream1) <- FALSE
rstream.sample(stream1,10)
