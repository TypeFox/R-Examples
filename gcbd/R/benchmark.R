
#hilbert <- function(n) {
#    i <- 1:n
#    1 / outer(i - 1, i, "+")
#}

getMatrix <- function(N) {
    a <- rnorm(N*N)
    #a <- hilbert(N)
    dim(a) <- c(N,N)
    invisible(gc())
    invisible(a)
}

## getMagmaMatrix <- function(N) {
##     a <- magma(rnorm(N*N), N, N)
##     #a <- magma( hilbert(N), N, N, gpu=TRUE)
##     invisible(gc())
##     invisible(a)
## }

matmultBenchmark <- function(N, n, trim=0.1) {
    a <- getMatrix(N)
    traw <- replicate(n, system.time(crossprod(a))[3])
    tmean <- mean(traw,trim=trim)
}

## matmultBenchmarkmagma <- function(N, n, trim=0.1) {
##     a <- getMagmaMatrix(N)
##     traw <- replicate(n, system.time(crossprod(a))[3])
##     tmean <- mean(traw,trim=trim)
## }

matmultBenchmarkgputools <- function(N, n, trim=0.1) {
    a <- getMatrix(N)
    traw <- replicate(n, system.time(gpuMatMult(a,a))[3])
    tmean <- mean(traw,trim=trim)
}

qrBenchmark <- function(N, n, trim=0.1) {
    a <- getMatrix(N)
    traw <- replicate(n, system.time(qr(a, LAPACK=TRUE))[3])
    tmean <- mean(traw,trim=trim)
}

## qrBenchmarkmagma <- function(N, n, trim=0.1) {
##     a <- getMagmaMatrix(N)
##     traw <- replicate(n, system.time(qr(a))[3])
##     tmean <- mean(traw,trim=trim)
## }

qrBenchmarkgputools <- function(N, n, trim=0.1) {
    a <- getMatrix(N)
    traw <- replicate(n, system.time(gpuQr(a))[3])
    tmean <- mean(traw,trim=trim)
}

svdBenchmark <- function(N, n, trim=0.1) {
    a <- getMatrix(N)
    traw <- replicate(n, system.time(svd(a))[3])
    tmean <- mean(traw,trim=trim)
}

#svdBenchmarkTEST <- function(N, n, trim=0.1) {
#    a <- getMatrix(N)
#    traw1 <- replicate(n, system.time(svd(a))); print(traw1)
#    traw <- replicate(n, system.time(svd(a))[3])
#    tmean <- mean(traw,trim=trim)
#}

## svdBenchmarkmagma <- function(N, n, trim=0.1) {
##     a <- getMagmaMatrix(N)
##     traw <- replicate(n, system.time(svd(a))[3])
##     tmean <- mean(traw,trim=trim)
## }

svdBenchmarkgputools <- function(N, n, trim=0.1) {
    a <- getMatrix(N)
    traw <- replicate(n, system.time(gpuSvd(a))[3])
    tmean <- mean(traw,trim=trim)
}

luBenchmark <- function(N, n, trim=0.1) {
    a <- getMatrix(N)
    traw <- replicate(n, system.time(lu(a))[3])
    tmean <- mean(traw,trim=trim)
}

## luBenchmarkmagma <- function(N, n, trim=0.1) {
##     a <- getMagmaMatrix(N)
##     traw <- replicate(n, system.time(lu(a))[3])
##     tmean <- mean(traw,trim=trim)
## }

luBenchmarkgputools <- function(N, n, trim=0.1) {
    NULL                                  # no LU in gputools
}
