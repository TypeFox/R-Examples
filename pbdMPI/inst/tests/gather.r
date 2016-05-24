suppressPackageStartupMessages(library(pbdMPI, quietly=TRUE))
init()

testval <- length(allgather(comm.rank()))
trueval <- comm.size()

if (comm.rank()==0) testval == trueval

finalize()
