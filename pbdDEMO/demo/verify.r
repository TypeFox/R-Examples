library(pbdDEMO, quietly=TRUE)
init.grid()

# set independent seeds using Rlecuyer
comm.set.seed(seed=1234, diff = TRUE)

# verify system solving at scale
verify.solve(nrows=1e3)

finalize()
