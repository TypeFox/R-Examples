library(pbdTEST)
settings(mpi=TRUE)

.BLDIM <- 2
comm.set.seed(seed=1234, diff=FALSE)


### --------------------------------------
module("Redistribute")

blacs_gridinit(ICTXT=3, NPROW=2, NPCOL=1L)

dx <- ddmatrix(1:100, 10, 10, bldim=c(6, 10), ICTXT=3L)

dy <- redistribute(dx, bldim=c(3,3), ICTXT=0)


test("Custom grid", {
  a <- as.matrix(dx)
  b <- as.matrix(dy)
})

collect()



finalize()

