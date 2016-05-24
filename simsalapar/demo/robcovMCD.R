require("simsalapar")

if(!require("robustbase"))
    stop("Needs the robustbase package; please install it, typically via\n",
	 "   install.packages(\"robustbase\")\n")
if(package_version(packageDescription("robustbase")$Version) < "0.9.8"
   || !any("raw.only" == names(formals(covMcd))))
    stop("Need a version of 'robustbase' where covMcd(*,raw.only=TRUE) works")

sessionInfo()# see it interactively
if(!exists("doExtras")) # <-- so you can set it *before*  demo(robcovMCD)
    doExtras <- simsalapar:::doExtras()
doExtras


## Here, we will redo part of the simulations done (and reported on) in
##
## Pison, G., Van Aelst, S., and Willems, G. (2002),
## Small Sample Corrections for LTS and MCD,
## _Metrika_ *55*, 111--123.
##
## ~/save/papers/robust-diverse/Pison_VanAelst_Willems.pdf

## build grid (for the *real* simulation
varList <- varlist(n.sim = list(type="N", expr = quote(N[sim]),
			       value = 1024),
		   ##		       ---- a multiple of block.size = 64 (below)
		   n     = list(type="grid", expr = quote(n),
				value = c(10, 15, 20, 30, 35, 40, 60, 100, 200)),
		   ## Note:  p < n is really "practically" required ...
		   p     = list(type="grid", expr = quote(p),
				value = c(1:5, 7, 10, 15, 25)),
		   alpha = list(type="grid", expr = quote(alpha),
				value = 1/2 + (0:8)/16)
		   )

dim(pGrid <- mkGrid(varList))
head(pGrid)
tail(pGrid)

.fullSim  <- exists(".fullSim")  && .fullSim   # default FALSE
.checking <- exists(".checking") && .checking  # default FALSE
if(!.fullSim && (interactive() || .checking)) {
    ## e.g., when run as a demo, need much smaller sizes
    message("Redimensioning the variable list in order to finish  ``in time'' ..")
    varList.full <- varList
    varList <- set.n.sim(varList, 512)
    varList$ n    $value <- local({n <- varList$n$value; n[n <= 60]})
    varList$ p    $value <- c(1:4, 7)
    varList$ alpha$value <- c(0.5, 0.75)
    print(pGrid <- mkGrid(varList))
}

ng <- get.nonGrids(varList)

## do1mcd() : to be applied to one line of the grid
do1mcd <- function(n, p, alpha)
{
    stopifnot(length(alpha) == 1, 0.5 <= alpha, alpha <= 1,
	      length(n) == 1, n == round(n), n > 0,
	      length(p) == 1, p == round(p), p > 0)
    X <- matrix(rnorm(n*p), n,p)
    ## Compute  Det(\hat\Sigma) ^ {1/p} :
    mcd <- covMcd(X, alpha=alpha, use.correction=FALSE, raw.only=TRUE)
    lD <- determinant(mcd$raw.cov)
    exp(as.vector(lD$modulus) / p) ## == Det(.)^{1/p}
}

##---- apply  do1mcd()  {*not* in parallel, for now}: ---------------------------

do1mcd(n=10, p=1, alpha=0.5)
do1mcd(n=30, p=2, alpha=0.5)

## a dummy small list, for testing -- this one with some p >= n:
vl.sm0 <- varlist(n.sim = list(expr = quote(N[sim]), value = 8),
		  n     = list(type="grid", expr= quote(n), value= 5*(2:4)),
		  ## NB: p >= n gives errors {but we can deal w/ them, after all!}
		  p     = list(type="grid", expr= quote(p), value = c(1:2,10,15)),
		  alpha = list(type="grid", expr= quote(alpha),value = c(.5, .75))
		  )
mkGrid(vl.sm0)# (as basic check) -> 24 rows

if(interactive())
    options(error = recover)
system.time(
r.sm0 <- doLapply(vl.sm0, seed = "seq", subjob.=subjob,
		  doOne=do1mcd, timer=mkTimer(gcFirst=FALSE))
)# 2.2 (lynne 2013)
## currently *with* errors :
str(e <- getArray(r.sm0,"error"))
print.table(e2 <- apply(e, 1:2, sum), zero=".")# error exactly  iff  p >= n
stopifnot(identical(which(e2 > 0), c(7L, 10:11)))

## a too small list; for testing (type = "N" is automatic for 'n.sim'):
vl.sml <- varlist(n.sim = list(expr = quote(N[sim]), value = 64),
		  n     = list(type="grid", expr = quote(n),
                               value = c(10, 15, 20:24, 30)),
		  p     = list(type="grid", expr = quote(p), value = c(1:3, 5)),
		  alpha = list(type="grid", expr = quote(alpha),
                               value = 1/2 + (0:2)/4)
		  )

smlFile <- "robcovMCD-sml-sim.rda"
## If the smlFile does not exist _or_ is older than 12 hours, (re)create it :
if(!file.exists(smlFile) ||
   Sys.time() - file.info(smlFile)[["ctime"]] > as.difftime(12, units="hours"))
{
    cat("recomputing ", smlFile, ": ") ; flush.console()
    stime.sml <- system.time(
        r.sml <- doLapply(vl.sml, seed = "seq",
                          subjob. = subjob, doOne=do1mcd, timer=mkTimer(gcFirst=FALSE)))
    print(stime.sml)# ~ 12 sec (lynne 2015-12); 23 sec (lynne 2014);  MH: 68s; ada-13: 45.5 sec; nb-mm3: 35.5 sec)
    save(stime.sml, r.sml, vl.sml, file = smlFile)
} else {
    cat("using", smlFile, "from" , format(file.info(smlFile)[["ctime"]]),":\n")
    load(smlFile)
}

(v.sml <- getArray(r.sml))
## Bug: these two plots are already contradicting; first plot is wrong: (p = 2: medians ~ 1.2)
mayplot(v.sml, vl.sml, row.vars = NULL, col.vars = "alpha", xvar = "n")
mayplot(v.sml, vl.sml, row.vars = "p",  col.vars = "alpha", xvar = "n", cex = .5)


if(! .checking) { ## now the "big" simulation:
    ##
print(system.time({
    res <- doClusterApply(varList, block.size = 64,  ## only 23 seonds
                 seed = "seq", #-> no 'seed' needed
                 initExpr = require("robustbase"), # <--!!{otherwise do1mcd() needs require(.): inefficient}
                 doOne=do1mcd, timer=mkTimer(FALSE))
}))
} else {
    cat("Using  **smaller**  example in the following:\n")
    res <- r.sml; varList <- vl.sml
}
## No warning/error anymore -- hooray

moreChecks <- FALSE
## NOTE: This works too --- about same speed on lynne 2015:
##  user  system elapsed
## 1.206   0.372  43.613
if(moreChecks)
print(system.time({
    reF <- doForeach(varList, block.size = 64,
                     seed = "seq", #-> no 'seed' needed
                     ## <<< 2014-12-04: takes 193 seconds on nb-mm3  ???
                     extraPkgs = "robustbase",
                     ## <--!!{otherwise do1mcd() needs require(.): inefficient}
                 doOne=do1mcd, timer=mkTimer(FALSE))
}))


if(!.checking) { ## save "all", including some platform info:
    print(sess.I <- sessionInfo())  # see it interactively
    node <- Sys.info()[c("nodename","release")]
    node.I <- c(node, cores = parallel::detectCores())
    sF <- sprintf("robcovMCD-sim_%s_%s.rda",
              node[1], format(Sys.Date()))
    cat(sprintf("Save file (in %s): %s\n", getwd(), sF))
    save(r.sml, vl.sml, res, varList, do1mcd, sess.I, node.I,
         file = sF)
}

## For interactive experiments
if(FALSE)
    load("robcovMCD-sim_ada-13_2014-11-27.rda")## takes ~ 2 sec !
if(FALSE)
    load("~/R/Pkgs/simsalapar.Rcheck-64b/tests/robcovMCD-sim_nb-mm3_2014-12-05.rda")

## convert array of nice simulated values
str(val <- getArray(res))
apply(val, 1:2, mean)## -- see boxplots below
apply(val, 1:2,   sd)## comparable

tools::assertError(## now gives an *intelligible* error message:
    ##   'row.vars' entry *not* in names(dimnames(x)): "meth"
mayplot(val, varList, pcol="tomato",
	row.vars = "meth", col.vars = "p", xvar = "n")
)

valm <- apply(val, 1:3, median)
dQ <- valm[,, if(.fullSim)"0.7500" else "0.75"] -
      valm[,, if(.fullSim)"0.5000" else "0.50"]
dQ  # alpha = 0.75 is consistently larger (medians!)
if(FALSE)## bug in mayplot() -- this gives an error {wrong switch() decision?}
mayplot(valm, varList, row.vars = NULL, col.vars = "p", xvar = "n", method="lines")

(med.dQ <- apply(dQ, 2, median))
## .checking == TRUE:
##         1         2         3         5
## 0.1776142 0.1378879 0.1333007 0.1365668
## .checking == FALSE -- .fullSim == FALSE
##         1         2         3         4         7
## 0.1278584 0.1244483 0.1235473 0.1275327 0.1108431
## .checking == FALSE -- .fullSim == TRUE
##         1         2         3         4         5         7  10 15 25
## 0.1389249 0.1179094 0.1134679 0.1154377 0.1176342 0.1022914  NA NA NA
if(.checking) {
    stopifnot(0.13 <= med.dQ, med.dQ <= 0.18)
} else if(!.fullSim) {
    stopifnot(0.11 <= med.dQ, med.dQ <= 0.13)
} else ## .fullSim
    stopifnot(0.10 <= med.dQ[1:6], med.dQ[1:6] <= 0.14)

## BUG in  mayplot():
mayplot(val, varList, row.vars = NULL, col.vars = "p", xvar = "n")
## looks "fine", but plotting wrong numbers: alpha = 0.75  medians > 1 is not ok !!!
round(valm[,, "0.75"], 2) ## where all are montonely *increasing* in 'n'
## which *is* clearly visible here [using new 'cex']:
mayplot(val, varList, pcol="tomato",
	row.vars = "alpha", col.vars = "p", xvar = "n", cex = .1, lwd=2)

### Using the dataframe and lattice:
dval <- array2df(val)
dvalm <- array2df(valm)

require(lattice)
## The medians only
xyplot(value ~ n | p + alpha, data=dvalm, type = c("p","smooth"))

## The full boxplots
bwplot(value ~ n | p + alpha, data=dval)

bwplot(value ~ n | p + alpha, data=dval, cex=0.75,
       xlab = "n", ylab = expression(Det(.) ^{1/p}),
       panel = function(...) { panel.grid(v = -1, h = -1); panel.bwplot(...) },
       par.settings = list(plot.symbol = list(cex = 0.25))) # smaller outlier symbols

str(times <- getArray(res, "time"))
table(round(times)) # of course quite platform dependent..

## check
stopifnot(val["20", "3", , 1:3] == res["20", "3", , 1:3]$value)

## Errors (trivial here, due to getArray()'s default FUN for comp="error"):
errs <- getArray(res, "error")
summary(errs)
n.err <- apply(errs, 1:2, sum) ## number of errors :
n.err                          ## ==> none at all for this example

## compute percentages of errors:
n.sim <- ng$n.sim
ftable(n.err/n.sim * 100) # all zero.. (as n.err ==== 0)
