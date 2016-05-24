require("simsalapar")

## build grid [type = "N" is automatic for 'n.sim']:
varList <- varlist(n.sim = list(expr = quote(N[sim]), value = 512),
                   n     = list(type="grid", expr = quote(n),
                                value = c(40, 100, 250)),
                   meth  = list(type="inner", expr = quote(italic(method)),
			        value = c("classical", "robust")))
## a version that has 'meth' as "grid" :
varL.1 <- varList
varL.1$meth$type <- "grid"

(pGrid <- mkGrid(varList))
ng <- get.nonGrids(varList)

## define doOne() (to be applied to one line of the grid)
## in a slightly modular version, and with extra argument 'DBG' :

##' Gaussian Mixture of .9 N(0,1) + .1 N(5, 2^2)
rMix <- function(n, pr2 = 0.1, mu2 = 5, sigma2 = 2) {
    stopifnot(is.numeric(n), n == round(n), length(n) == 1,
              0 <= pr2, pr2 <= 1, sigma2 >= 0, is.numeric(mu2))
    n2 <- rbinom(1, size=n, prob = pr2)
    sample(c(rnorm(n - n2), mu2 + sigma2*rnorm(n2)))
}

## meth="grid" (==> no nonGrids) -- with 'varL.1' :
do.1 <- function(n, meth, DBG=FALSE) {
    stopifnot(length(meth) == 1)
    y <- rMix(n)
    ## Now, simply
    ##	 if(meth == "robust") mean(y, trim = .25) else mean(y)
    ## but in a way to use  doCallWE(.) :
    M <- if(identical(meth, "robust")) function(u) mean(u, trim=.25) else mean
    r <- M(y)
    if(DBG)
        cat(sprintf("do1(x): n=%4d, meth=%9s : -> value %12g\n",
                    n, meth, r))
    r
}
codetools::checkUsage(do.1)
## codetools::findGlobals(do.1, FALSE)$variables ## "empty"

## faster: meth="inner" , with 'varList' :
do.one <- function(n, meth, DBG=FALSE) {
    stopifnot(is.character(meth))
    y <- rMix(n)
    M <- function(u, method) if(method == "robust") mean(u, trim=.25) else mean(u)
    ## On *this* data, use all methods:
    r <- vapply(meth, M, numeric(1), u=y)
    if(DBG)
	cat(sprintf("doOne(x): n=%4d : -> values %8.4f %8.4f\n",
		    n, r[1], r[2]))
    r
}
codetools::checkUsage(do.one)

## Quickly demonstrate what's going on.
## Note this is *fast* as we do not do any error catching !
system.time({
set.seed(1); rr <- replicate(1000, mean(rMix(100), trim=.25))
set.seed(1); rc <- replicate(1000, mean(rMix(100)))
}) # 0.2 seconds

boxplot(cbind(mean = rc, 'mean(trim = 0.25)' = rr),
        notch = TRUE); abline(h = 0, lty=3, col="gray20")

##---- test do.one() :
do.one(n=50, meth=c("robust","classical"))

##---- apply  do.one()  {*not* in parallel, for now}: ---------------------------

## Demonstrating the DBG extra argument:
vl12 <- set.n.sim(varList, 12)
system.time(
r12 <- doLapply(vl12, DBG=TRUE,
               doOne=do.one, timer=mkTimer(gcFirst=FALSE))
)
(v12 <- getArray(r12))
mayplot(v12, vl12, pcol="tomato", row.vars = NULL, col.vars = "meth", xvar = "n",
        notch = TRUE)

## If system.time() is called via its default, gcFirst=TRUE,
## garbage collection uses up 99% of the time (in this fast example):
system.time(
r12g <- doLapply(vl12, monitor=TRUE, # <- generally provided alternative to 'DBG'
                 doOne=do.one, timer=mkTimer(gcFirst=TRUE))
)
## (compare with the much shorter time for 'r12' above)
stopifnot(doRes.equal(r12, r12g))

## Compare the meth = "grid" with meth = "inner":
v.12 <- set.n.sim(varL.1, 12)
system.time(
g12 <- doLapply(v.12, monitor=TRUE, doOne=do.1, timer=mkTimer(gcFirst=TRUE))
)
stopifnot(all.equal(      getArray(g12),
                    aperm(getArray(r12), c(2:1,3)),
		    check.attributes=FALSE))## <- buglet in R: aperm(a) loses names(dim(a))


(tL <- system.time({
res <- doLapply(varList, doOne=do.one, timer=mkTimer(FALSE))
}))
## quite fast (*because* turned off gc in timer!)

## convert to an array of lists
## str(res <- mkAL(r0, vList=varList), max=0)
str(res, max=0)


## convert array of nice simulated values
str(val <- getArray(res))
apply(val, 1:2, mean)## -- see boxplots below
apply(val, 1:2,   sd)## comparable

mayplot(val, varList, pcol="tomato",
        ## TODO MM: these should be default; no "inner" -> take first "grid"
        row.vars = NULL, col.vars = "meth", xvar = "n")
if(FALSE)## fails (clearly), but one could think it should work (1 x 1 panel; inner var 'v')
mayplot(val, varList, pcol="tomato",
        row.vars = NULL, col.vars = NULL, xvar = "n")

str(times <- getArray(res, "time"))
table(round(times)) # most are '0'; rest very platform dependent..
## still, get average time in milli-seconds:
round(1000*apply(times,1:2, mean))

## check
stopifnot(identical(val["robust", "100", 1:3],
                    vapply(res["100", 1:3], `[[`, numeric(2),
                           "value")["robust",]))

## Errors (trivial here, due to getArray()'s default FUN for comp="error"):
errs <- getArray(res, "error")
summary(errs)
n.err <- apply(errs, 1:2, sum) ## number of errors :
n.err                          ## ==> none at all for this example

## compute percentages of errors:
n.sim <- ng$n.sim
ftable(n.err/n.sim * 100)

## Now in parallel: -----------------------------------

## due to 'R CMD check --as-cran' allowing only <= 2 cores
(nc <- simsalapar:::nCores4test()) # full cores for .parallel.chk.users or via env.var

do1 <- FALSE  ## too slow: 24-30 seconds (on lynne)
if(do1)
(tM1 <- system.time({
resC <- doMclapply(varList, cores=nc, block.size = 16,
		   doOne=do.one, timer=mkTimer(FALSE))
}))

## This time using blocks -- using a "lapply loop" :
bl.sizes <- c(8,16,32,64,128,256)
bl.sizes <- setNames(bl.sizes, paste("bs", bl.sizes, sep="="))
rrr <- lapply(bl.sizes, function(bl.size) {
    cat("block.size = ", bl.size,"\n")
    st <- system.time(
        r <- doMclapply(varList, cores=nc, block.size = bl.size,
                        doOne=do.one, timer=mkTimer(FALSE)) )
    list(res = r, sysTime = st)
})
object.size(rrr)## 13.4 MB {now w/o .Random.seed's}

lrr <- lapply(rrr, `[[`, 1)

## Check equality of the important part of the results:
if(do1) stopifnot(doRes.equal(res, resC))
for(j in seq_along(bl.sizes)) {
    cat(format(names(bl.sizes))[j],": ",
        if(isTRUE(ee <- doRes.equal(res, lrr[[j]]))) "Ok"
        else paste("**:", ee), "\n", sep="")
}


## Which block size is fastest?   Look at "elapsed":
rbind(tL, if(do1) tM1, t(sapply(rrr, `[[`, "sysTime")))
## (the optimum depends quite a bit on your platform, #{cores}, ...

## Martin's nb-mm3; Lenovo X201; Linux 3.8.0-25-generic #37-Ubuntu SMP x86_64:
##        user.self sys.self elapsed user.child sys.child
## tL         1.676    0.024   1.705      0.000     0.000
## bs=8       0.236    0.988   2.196      3.184     2.948
## bs=16      0.136    0.504   1.701      3.072     2.304
## bs=32      0.064    0.264   1.368      3.024     1.656
## bs=64      0.048    0.136   1.295      3.028     1.400
## bs=128     0.020    0.080   1.124      2.768     0.868
## bs=256     0.016    0.044   0.887      2.496     0.536

## Martin's lynne; AMD Phenom II X4 925; Linux 3.8.13-100.fc17.x86_64:
##        user.self sys.self elapsed user.child sys.child
## tL         0.839    0.003   0.850      0.000     0.000
## bs=8       0.094    1.042   1.240      1.054     1.928
## bs=16      0.074    0.538   0.834      0.951     1.363
## bs=32      0.061    0.272   0.622      0.950     0.829
## bs=64      0.050    0.144   0.561      0.998     0.631
## bs=128     0.045    0.075   0.492      0.824     0.537
## bs=256     0.047    0.031   0.494      0.770     0.337


c(Sys.info()[c(4:5,1:2)], cores = parallel::detectCores())
