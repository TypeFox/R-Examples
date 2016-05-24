library(rexpokit, quietly=TRUE)
library(pbdDMAT, quietly=TRUE)
init.grid()


.pbd_env$SPMD.CT$print.quiet <- TRUE


n <- 6e2

if (comm.rank()==0){
  x <- matrix(rnorm(n*n), n, n)
} else {
  x <- NULL
}

dx <- as.ddmatrix(x)

if (comm.rank()==0){
    ts <- system.time(a <- expm(x))
    te <- system.time(d <- rexpokit::expokit_dgpadm_Qmat(x, t=1))
}

barrier()

tp <- system.time(db <- expm(dx))
b <- as.matrix(db)


comm.cat("\n\n")
comm.cat(sprintf("rexpokit: %f\nserial:   %f\nparallel: %f\n\n", te[3], ts[3], tp[3]))
#
comm.cat(sprintf("expokit==serial?   "))
comm.cat(all.equal(a, d))
comm.cat("\n")
#
comm.cat(sprintf("serial==parallel?  "))
comm.cat(all.equal(a, b))
comm.cat("\n")

finalize()




