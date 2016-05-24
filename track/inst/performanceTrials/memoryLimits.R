# process 40 16Mb objects (about 640Mb) with memory limited to 128Mb
library(trackObjs)
source("trackObjs/inst/performanceTrials/funsForTesting.R")
options(width=100)
sessionInfo()
mem.limits(nsize=2^20, vsize=2^27) # vsize=128Mb
mem.limits()
unlink("tmp1", recursive=TRUE)
track.start("tmp1", options=list(cache=FALSE))
totSize <- 0
for (i in 1:40) {
    x <- createTestObj(1, 14) # scale=14 is a 16Mb object
    n <- randObjName(3, simple=TRUE, max.len=7)
    track(list=n)
    assign(n, value=x)
    totSize <- totSize + object.size(x)
    cat("Done obj", i, "totSize=", totSize, "\n")
}
track.summary()
track.stop()
unlink("tmp1", recursive=TRUE)
