himed <- function(x) { n2 <- 1 + length(x) %/% 2; sort(x, partial = n2)[n2] }

## Tolerance  2e-7 {13 * higher than default  1.49e-8 }
is.all.equal <- function(x,y, tol = 2e-7) {
    is.logical(r <- all.equal(x,y, tolerance = tol)) && r }


library(robustbase)

options(digits = 7)# single precision!
set.seed(15)

cat("  n |   range(x)   | wgt.Himed\n",
    "------------------------------\n",sep="")
for(i in 1:100) {
    n <- rpois(1, lam = 10)
    cat(formatC(n,wid=3)," ")
    x <- round(rnorm(n),3)
    iw <- 1 + rpois(n, lam = 2)
    him   <- himed(rep(x, iw)) ## == naive R solution
    whim <- wgt.himedian (x, iw)
    if(!is.all.equal(whim, him))
        cat("whim != him:    ", whim, "!=", him,"\n")
    cat(formatC(range(x), wid = 6, flag="-"), "",
        formatC(whim,     wid = 6, flag="+"), "\n")
}

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
