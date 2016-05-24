stopifnot(require("Rmpfr"))

stopifnot(chooseMpfr(1:10, 0) == 1,# failed earlier
	  chooseMpfr(20, 0:20) == choose(20, 0:20),
	  chooseMpfr(19, 0:20) == choose(19, 0:20),
	  chooseMpfr	(30, 4:30) * (-1)^(4:30) ==
	  chooseMpfr.all(30, k0=4, alternating=TRUE)
          )

cat('Time elapsed: ', proc.time(),'\n') # "stats"

## sumBinomMpfr() ... had embarrasing bug for a while
sBn <- Rmpfr:::sumBinomMpfr.v1
stopifnot(
    all.equal(         sBn(10, sqrt),
              sumBinomMpfr(10, sqrt), tol=1e-77) ,
    all.equal(         sBn(10, log, n0=1, alternating=FALSE),
              sumBinomMpfr(10, log, n0=1, alternating=FALSE), tol=1e-77)
    )

fBin <- function(k) x^k * (1-x)^(n-k)
## \sum_{k=0}^n  (n \\ k) x^k (1-x)^{n-k} == sum(dbinom(0:n, n, prob=x)) == 1 :
for(x in runif(50)) {
    n <- 1 + rpois(1, lambda=10)
    cat(".")
    stopifnot(all.equal(1, sumBinomMpfr(n, fBin, alternating=FALSE),
                        tol = 1e-15))
};cat("\n")


cat('Time elapsed: ', proc.time(),'\n') # "stats"

if(!interactive()) warnings()
