suppressPackageStartupMessages(library(pbdMPI, quietly=TRUE))
init()

cs <- comm.size()

testval <- allreduce(comm.rank())
trueval <- cs*(cs-1)/2

test <- comm.all(testval == trueval)
comm.print(test)

finalize()
