
library(robustbase)
## testing functions:
source(system.file("xtraR/ex-funs.R", package = "robustbase"))

x  <- c(0.26, 0.161, 1.33, -0.925, 0.199, -1.476, 0.489)
iw <- c(5, 4, 4, 1, 5, 1, 5)


stopifnot(0.26 == (himR <-  weighted.median(rep(x,iw))),
          himR == wgt.himedian(x, iw), ## (once gave infinite loop)
          himR == wgt.himedian(x, as.integer(iw)))


## same result, but *different  wweigted.median() debug output!
##-- even when having EXACT data (& exact differences!)
all.equal(Qn(c(2:1,4:3)), 1.1376128)

###--- another inifinite loop {solved}:

(z4 <- round(qnorm(ppoints(4)), 2))

## both the same (also wweigted.median debug output)
(all.equal(weighted.median(z4, 4:1),
           print(wgt.himedian (z4, 4:1))))# 3.97e-8

(all.equal(weighted.median(z4, c(4,2,3,17)),
           print(wgt.himedian (z4, c(4,2,3,17)))))# 4.54e-8

Sn (z4)## = 0.8533053
Sn (z4, const = 1)# = 0.75

##-- now
Qn (z4)# --> gave (another) infinite loop
       ##--> now "works" after (float) rounding of differences!
       ##--- DIFFERENT whimed() output!
stopifnot(all.equal(Qn(z4, const = 1),
                    print(Qn0R(z4))))

## yet another problem:
Sn0R(c(1.1, -0.93, -0.11, -0.74))# 0.82
Sn  (c(1.1, -0.93, -0.11, -0.74))# 0.9329471
## gave segmentation fault at Sat Mar 16 23:54:30 2002
## not anymore but 0.9329471

### Check validity of basic algorithm few times
set.seed(471)
for(sim in 1:100) { # had '500'
    cat(".")
    x <- rnorm(rpois(1, lam=80))# not too large the *n0R() use time!
    ##--> Sn0R() "fails" for  odd n
    stopifnot(all.equal(Sn(x, const = 1), Sn0R(x)),
              all.equal(Qn(x, const = 1), Qn0R(x), tolerance = 7e-8))
    x <- round(x,2)
    stopifnot(all.equal(Sn(x, const = 1), Sn0R(x)),
              all.equal(Qn(x, const = 1), Qn0R(x), tolerance = 7e-8))
    if(sim %% 50 == 0) cat(sim, "\n")
}

###---- Last series of problems: when  n^2 > max.integer:

## Large x with 1% outliers
N <- 1e5
n.o <- round(0.01 * N)
nSim <- 24## interesting
nSim <- 4 ## for package testing
estim.lst <- c("mad", "Sn", "Qn")
Res <- array(NA, dim = c(nSim, length(estim.lst), 1 + 2),
             dimnames= list(NULL,estim.lst, c("Tx","cpu1", "cpu3")))
set.seed(101)
for(i in 1:nSim) {
    x <- sample(c(rnorm(N), 10*min(1, abs(rt(1, df=2))) + rnorm(n.o)))
    cat(i)
    for(S in estim.lst) {
        cpu <- system.time(Tx <- get(S)(x))[1:3]
        Res[i, S,] <- c(Tx, cpu[c(1,3)])
    }
    cat(" ")
};  cat("\n")

options(digits = 5)

(Tx <- Res[,, "Tx"])
stopifnot(abs(range(Tx - 1)) < 0.03)

q()

### -- Rest: rather for demo -- keep here for reference

apply(Res, c(2,3), mean)

## Variation: robust or not:
1000* apply(Tx, 2, sd)#-> Qn < Sn < mad
1000* apply(Tx, 2, Qn)#-> Qn > Sn > mad

if(dev.interactive(orNone=TRUE)) {
    boxplot(Tx, main = sprintf("n=%d  x N(0,1) + %d (1%%) outliers to the right",
                N,n.o))
    abline(h = 1, lty = 3, lwd = 2, col = "gray")
}

if(interactive()) { ## i.e. not when package testing ..

N <- 500
set.seed(101)
str(iw <- 1L+ as.integer(rpois(N, 1))); str(w <- as.double(iw))

cr <- ci <- numeric(50)
for(nn in seq_along(ci)) {
    x <- round(rnorm(N),1)
    cat(".")
    cr[nn] <- system.time(for(i in 1:1000) rr <- wgt.himedian(x,  w))[1]
    ci[nn] <- system.time(for(i in 1:1000) ri <- wgt.himedian(x, iw))[1]
    stopifnot(rr == ri)
};cat("\n")


## Or rather (as correctly) a "paired" comparsion:
boxplot(cr - ci, notch=TRUE)
## rather
t.test(     cr,      ci, paired = TRUE) ##-> P-value of 0.0219
t.test(log(cr), log(ci), paired = TRUE) ##-> P-value of 0.0088
wilcox.test(cr, ci, paired = TRUE)      ##-> P-value of 2.23e-5 (!!)

}
