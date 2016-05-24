library(expm)

source(system.file("test-tools.R", package= "expm"), keep.source=FALSE)# relErr()

set.seed(101)
for(n in c(1:5, 10:11, if(doExtras) 100:101 else 25)) {
    cat("n = ",n,"\n-----\n")
    for(i in seq_len(if(doExtras)10 else 3)) {
        A <- matrix(round(10*rnorm(n^2))/4, n,n)
        E <- matrix(rnorm(n^2, sd = 1e-3),  n,n)
        F1 <- expmFrechet(A, E)
        F2 <- expmFrechet(A, E, "block")
        if(i == 1 && n < 9) print(F1)
        stopifnot(all.equal(F1, F2, tol = 6e-15 * n))
        cat(sprintf("%5.2f ", relErr(F1 $ L, F2 $ L) * 2^52))
    }
    cat(" * eps_C \n")
}

cat('Time elapsed: ', proc.time(), '\n') # for "statistical reasons"
