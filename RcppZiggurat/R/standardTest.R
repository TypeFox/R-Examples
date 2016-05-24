## This follows the approach in (R)DieHarder which does
##   take N draws from the N(0,1) we test here
##   convert into U(0,1) by using the inverse of the normal
##   repeat M times
##   and for large enough N, then the sum of all N draws goes to
##       mean   --> N/2
##       stddev --> sqrt(N/12)
##   which is known as the Irwin-Hall distribution
##   then for each of these M values use the inverse of normal to obtain a p-value
##   that p value should be uniformly distributed across these M draws
##   so use Kuiper's K/S test variant to test for uniform U(0,1)

standardTest <- function(N=1e5,      	# individual draws
                         M=1e2,  	# repeats
                         seed=123456789,
                         generators=c("Ziggurat", "MT", "LZLLV",
                                      "GSL", "QL", "Gretl"),
                         showplot=interactive()) {

    res <- mclapply(generators, FUN=function(g, seed) {
        res <- ziggtest(M, N, g, seed)
        v <- pnorm(res, mean=N/2, sd=sqrt(N/12))
    }, seed, mc.cores=getOption("mc.cores", 2L))
    names(res) <- generators
    res <- as.data.frame(res)

    attr(res, "testtype")  <- "Standard"
    attr(res, "draws")     <- N
    attr(res, "repeats")   <- M
    attr(res, "seed")      <- seed
    attr(res, "created")   <- format(Sys.time())
    attr(res, "version")   <- packageVersion("RcppZiggurat")

    if (showplot) {
        plotAll(res)
    }

    invisible(res)
}
