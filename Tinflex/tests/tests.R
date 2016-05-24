#############################################################################
##                                                                         ## 
## Tests for Tinflex.                                                      ## 
##                                                                         ## 
#############################################################################

## Load libraries -----------------------------------------------------------

require(Tinflex)

## Constants ----------------------------------------------------------------

## Seed for running tests.
SEED <- 123456
set.seed(SEED)

## Sample size for goodness-of-fit (GoF) tests.
## (if n.GoF, then is >= 1e7 then it is rounded to a multiple of 1e7.)
n.GoF <- 1e5

## Sample size for creating histogram.
## (if <= 0, then no histgram is plotted.)
n.hist <- 1e4

## Sample size for comparing R and C version of sampling routine.
## (if <= 0, then test is skipped.)
n.comp <- 100

## Requested maximal ratio rho = A.hat / A.squeeze.
rho <- 1.1

## Lower bounded for accepted p-values.
pval.threshold <- 1e-4

## Print numbers very accurately.
options(digits = 15)

## Load auxiliary libraries -------------------------------------------------

if (  
    ## Routines for testing non-uniform random variate generators.
    require(rvgtest) &&
    ## Function for approximate quantile function.
    require(Runuran) ) {

  ## We have all required libraries.
  have.UNURAN <- TRUE

} else {
  warning("Packages 'Runuran' and 'rvgtest' needed for testing!")
  have.UNURAN <- FALSE
}
 
## Global variables ---------------------------------------------------------

## Start timer.
time.start <- proc.time()

## Vector for storing p-values of tests.
pvals <- numeric(0)

## Number of compare tests that fail.
comp.fail <- 0

## Auxiliary functions ------------------------------------------------------

## Plot histogram and run GoF tests (only if 'Runuran' is installed). .......
run.test <- function(gen, cT, lf, lb=-Inf, ub=Inf, plot=FALSE) {
  ## Print generator in debugging mod.e
  print(gen, debug=TRUE)

  ## Plot density, hat and squeeze.
  if (plot) {
    plot(gen, from=max(lb,-3), to=min(ub,3), is.trans=FALSE, main=paste("c =",cT))
    plot(gen, from=max(lb,-3), to=min(ub,3), is.trans=TRUE,  main=paste("c =",cT))
  }
  
  ## Compare R and C version of sampling routine.
  if (n.comp > 0) {
    set.seed(SEED)
    x.R <- Tinflex:::Tinflex.sample.R(gen,n=n.comp)
    set.seed(SEED)
    x.C <- Tinflex.sample(gen,n=n.comp)
    if (!identical(x.R,x.C)) {
      comp.fail <<- comp.fail + 1
      d <- x.R - x.C
      cat(paste("Warning:!\n",
                "R and C versions differ!\n",
                "\tsample size =",n.comp,"\n",
                "\t  different =",length(d[d!=0]),"\n\n"))
      cat("output of R version:\n")
      print(x.R[d!=0])
      cat("difference to C version:\n")
      print(d[d!=0])
      ## Remark: x.R and x.C may differ due to different round-off errors
      ## in the C and R version, resp.
      ## Thus we only print a warning.
    }
  }
  
  ## Create a histgram and run GoF test (when Runuran is available).
  if (isTRUE(have.UNURAN)) {
    ## Set seed for tests.
    set.seed(SEED)

    ## Create frequency table (using packages 'Runuran' and 'rvgtest')
    if (n.GoF < 1e7) {
      m <- 1
      n <- n.GoF
    } else {
      m <- round(n.GoF / 1e7)
      n <- 1e7
    }      
    ug <- pinv.new(pdf=lf, lb=lb, ub=ub, islog=TRUE, center=0, uresolution=1.e-14)
    tb <- rvgt.ftable(n=n, rep=m,
                      rdist=function(n){Tinflex.sample(gen,n=n)},
                      qdist=function(u){uq(ug,u)})

    ## Plot histgram of random sample.
    if (plot && n.hist > 0)
      hist(Tinflex.sample(gen,n=n.hist), breaks=101, main=paste("c =",cT))

    ## Plot frequency table,
    if (plot) {
      plot(tb, main=paste("c =",cT))
    }

    ## Run GoF test.
    gof <- rvgt.chisq(tb)
    print(gof)
    if(! isTRUE(gof$pval[m] >= pval.threshold))
      warning("p-value too small. GoF test failed!")

    ## Store p-value.
    pvals <<- append(pvals, gof$pval[m])
  }
  else {
    ## Package 'Runuran' or 'rvgtest' is not installed.
    ## Thus we only plot a histgram.
    if (plot && n.hist > 0)
      hist(Tinflex.sample(gen,n=n.hist), breaks=101, main=paste("c =",cT))
  }
}

## Test whether there is an error. ..........................................
is.error <- function (expr) { is(try(expr), "try-error") }


#############################################################################
## Check API
#############################################################################

lf <- function(x) { -x^4 + 5*x^2 - 4 }  ## = (1 - x^2) * (x^2 - 4)
dlf <- function(x) { 10*x - 4*x^3 }
d2lf <- function(x) { 10 - 12*x^2 }

if (! (is.error( Tinflex.setup() ) &&
       is.error( Tinflex.setup(lpdf=1) )) )
    stop("invalid 'lpdf' not detected")

if (! (is.error( Tinflex.setup(lpdf=lf) ) &&
       is.error( Tinflex.setup(lpdf=lf, dlpdf=1) )) )
    stop("invalid 'dlpdf' not detected")

if (! (is.error( Tinflex.setup(lpdf=lf, dlpdf=dlf) ) &&
       is.error( Tinflex.setup(lpdf=lf, dlpdf=dlf, d2lpdf=1) )) )
    stop("invalid 'd2lpdf' not detected")

if (! (is.error( Tinflex.setup(lpdf=lf, dlpdf=dlf, d2lpdf=d2lf) ) &&
       is.error( Tinflex.setup(lpdf=lf, dlpdf=dlf, d2lpdf=d2lf, ib=1) )) )
    stop("invalid 'ib' not detected")

if (! (is.error( Tinflex.setup(lpdf=lf, dlpdf=dlf, d2lpdf=d2lf, ib=c(-Inf,0,Inf), cT="a") ) &&
       is.error( Tinflex.setup(lpdf=lf, dlpdf=dlf, d2lpdf=d2lf, ib=c(-Inf,0,Inf), cT=c(-0.5,0,-0.5)) )) )
    stop("invalid 'cT' not detected")

if (! (is.error( Tinflex.setup(lpdf=lf, dlpdf=dlf, d2lpdf=d2lf, ib=c(-Inf,0,1,Inf), cT=c(-2,0,-0.5)) ) &&
       is.error( Tinflex.setup(lpdf=lf, dlpdf=dlf, d2lpdf=d2lf, ib=c(-Inf,0,1,Inf), cT=c(0,-0.5,-2)) )) )
    stop("invalid 'cT' not detected")


#############################################################################
## Distribution 1
#############################################################################

lf <- function(x) { -x^4 + 5*x^2 - 4 }  ## = (1 - x^2) * (x^2 - 4)
dlf <- function(x) { 10*x - 4*x^3 }
d2lf <- function(x) { 10 - 12*x^2 }

## extrema: -1.581, 0.0, 1.581

## c = 1.5 ------------------------------------------------------------------
## inflection points: -1.7620, -1.4012, 1.4012, 1.7620
cT <- 1.5

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-3,-1.5,0,1.5,3), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3,ub=3)

## c = 1 --------------------------------------------------------------------
## inflection points: -1.8018, -1.3627, 1.3627, 1.8018
cT <- 1

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-3,-1.5,0,1.5,3), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3,ub=3)

## c = 0.5 ------------------------------------------------------------------
## inflection points: -1.8901, -1.2809, 1.2809, 1.8901
cT <- 0.5

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-3,-1.5,0,1.5,3), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3,ub=3)

## c = 0.1 ------------------------------------------------------------------
## inflection points: -2.2361, -1.0574, 1.0574, -2.2361
cT <- 0.1

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-3,-1.5,0,1.5,3), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3,ub=3)

## c = 0 --------------------------------------------------------------------
## inflection points: -0.9129, 0.9129
cT <- 0

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-Inf,-2.1,-1.05,0.1,1.2,2,Inf), cT=cT, rho=rho)
run.test(gen,cT,lf)

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-Inf,-1,0,-1,Inf), cT=cT, rho=rho)
run.test(gen,cT,lf)

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-2,0,1.5), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-2,ub=1.5)

## c = -0.2 -----------------------------------------------------------------
## inflection points: -0.4264, 0.4264
cT <- -0.2

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-Inf,-2.1,-1.05,0.1,1.2,2,Inf), cT=cT, rho=rho)
run.test(gen,cT,lf)

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-Inf,-1,0,1,Inf), cT=cT, rho=rho)
run.test(gen,cT,lf)

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-2,0,1.5), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-2,ub=1.5)

## c = -0.5 -----------------------------------------------------------------
## inflection points: -0.4264, 0.4264
cT <- -0.5

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-Inf,-2.1,-1.05,0.1,1.2,2,Inf), cT=cT, rho=rho)
run.test(gen,cT,lf)

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-Inf,-1,0,1,Inf), cT=cT, rho=rho)
run.test(gen,cT,lf)

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-2,0,1.5), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-2,ub=1.5)

## c = -0.9 -----------------------------------------------------------------
## inflection points: -0.4264, 0.4264
cT <- -0.9

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-Inf,-2.1,-1.05,0.1,1.2,2,Inf), cT=cT, rho=rho)
run.test(gen,cT,lf)

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-Inf,-1,0,1,Inf), cT=cT, rho=rho)
run.test(gen,cT,lf)

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-2,0,1.5), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-2,ub=1.5)

## c = -1 -------------------------------------------------------------------
## inflection points: -0.3094, -0.3094
cT <- -1

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-3,-2.1,-1.05,0.1,1.2,3), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3,ub=3)

## c = -1.5 -----------------------------------------------------------------
## inflection points: -0.2546, 0.2546
cT <- -1.5

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-3,-2.1,-1.05,0.1,1.2,3), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3,ub=3)

## c = -2 -------------------------------------------------------------------
## inflection points: -0.2213, 0.2213
cT <- -2

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-3,-2.1,-1.05,0.1,1.2,3), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3,ub=3)


#############################################################################
## Distribution 2
## Test construction points at and near extrema.
#############################################################################

lf <- function(x) { -2*x^4 + 4*x^2 } 
dlf <- function(x) { -8*x^3 + 8*x }
d2lf <- function(x) { -24*x^2+8 }

## extrema: -1, 0, 1

## special boundary points:
##   points where the slope of the tangent is identical 0.
ivb.1 <- c(-Inf, -2, -1, 0, 1, 2, Inf) 
##   points where the slope of the tangent is almost 0
##   but not identical 0 (and thus might cause numerical errors).
ivb.2 <- c(-Inf, -2, -1+2^(-52), 1e-20, 1-2^(-53), 2, Inf) 

ivb.3 <- c(-3, -1, 0, 1, 3) 
ivb.4 <- c(-3, -1+2^(-52), 1e-20, 1-2^(-53), 3) 

## c = 2 --------------------------------------------------------------------
## inflection points: -1.1734, -0.8300, 0.8300, 1.1734
cT <- 2

## slope = 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.3, cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3, ub=3)

## slope ~ 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.4, cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3, ub=3)
rm(gen)

## c = 1.5 ------------------------------------------------------------------
## inflection points: -1.1993, -0.8067, 0.8067, 1.1993
cT <- 1.5

## slope = 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.3, cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3, ub=3)

## slope ~ 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.4, cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3, ub=3)

## c = 1 --------------------------------------------------------------------
## inflection points: -1.2418, -0.7709, -0.7709, 1.2418
cT <- 1

## slope = 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.3, cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3, ub=3)

## slope ~ 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.4, cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3, ub=3)

## c = 0.5 ------------------------------------------------------------------
## inflection points: -1.3345, -0.7071, 0.7071, 1.3345
cT <- 0.5

## slope = 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.3, cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3, ub=3)

## slope ~ 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.4, cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3, ub=3)

## c = 0 --------------------------------------------------------------------
## inflection points: -0.5774, 0.5774, 
cT <- 0

## slope = 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.1, cT=cT, rho=rho)
run.test(gen,cT,lf)

## slope ~ 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.2, cT=cT, rho=rho)
run.test(gen,cT,lf)

## c = -0.2 -----------------------------------------------------------------
## inflection points: -0.5076, 0.5076
cT <- -0.2

## slope = 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.1, cT=cT, rho=rho)
run.test(gen,cT,lf)

## slope ~ 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.2, cT=cT, rho=rho)
run.test(gen,cT,lf)

## c = -0.5 -----------------------------------------------------------------
## inflection points: -0.4180, 0.4180
cT <- -0.5

## slope = 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.1, cT=cT, rho=rho)
run.test(gen,cT,lf)

## slope ~ 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.2, cT=cT, rho=rho)
run.test(gen,cT,lf)

## c = -0.9 -----------------------------------------------------------------
## inflection points: -0.3404, 0.3404
cT <- -0.9

## slope = 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.1, cT=cT, rho=rho)
run.test(gen,cT,lf)

## slope ~ 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.2, cT=cT, rho=rho)
run.test(gen,cT,lf)

## c = -1 -------------------------------------------------------------------
## inflection points: -0.3264, 0.3264
cT <- -1

## slope = 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.3, cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3,ub=3)

## slope ~ 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.4, cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3,ub=3)

## c = -1.1 ------------------------------------------------------------------
## inflection points: -0.3139, 0.3139
cT <- -1.1

## slope = 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.3, cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3,ub=3)

## slope ~ 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.4, cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3,ub=3)

## c = -2 --------------------------------------------------------------------
## inflection points: -0.2412, 0.2412
cT <- -2

## slope = 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.3, cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3,ub=3)

## slope ~ 0
gen <- Tinflex.setup(lf, dlf, d2lf, ib=ivb.4, cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3,ub=3)

rm(ivb.1, ivb.2, ivb.3, ivb.4)


#############################################################################
## Distribution 3
## Test density that vanishes at boundaries of bounded domain.
#############################################################################

lf <- function(x) { log(1-x^4) }
dlf <- function(x) { -4*x^3/(1-x^4) }
## extrema: 0
d2lf <- function(x) { -(4*x^6+12*x^2)/(x^8-2*x^4+1) } 

## domain: lb=-1, ub=1

## c = 2   ------------------------------------------------------------------
## inflection points: -0.8091, [ 0 ], 0.8091
cT <- 2

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-1,-0.9,-0.5,0.5,0.9,1), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-1,ub=1)

## c = 1.5 ------------------------------------------------------------------
## inflection points: -0.8801, [ 0 ], 0.8801
cT <- 1.5

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-1,-0.9,-0.5,0.5,0.9,1), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-1,ub=1)

## c = 1 --------------------------------------------------------------------
## inflection points: [ 0 ]
cT <- 1

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-1,-0.5,0,0.5,1), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-1,ub=1)

## c = 0.5 ------------------------------------------------------------------
cT <- 0.5

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-1,-0.5,0,0.5,1), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-1,ub=1)

## c = 0.1 ------------------------------------------------------------------
cT <- 0.1

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-1,-0.5,0,0.5,1), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-1,ub=1)

## c = 0 --------------------------------------------------------------------
cT <- 0

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-1,-0.5,0,0.5,1), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-1,ub=1)

## c = -0.2 -----------------------------------------------------------------
cT <- -0.2

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-1,-0.5,0,0.5,1), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-1,ub=1)

## c = -0.5 -----------------------------------------------------------------
cT <- -0.5

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-1,-0.5,0,0.5,1), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-1,ub=1)

## c = -0.9 -----------------------------------------------------------------
cT <- -0.9

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-1,-0.5,0,0.5,1), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-1,ub=1)

## c = -1 -------------------------------------------------------------------
cT <- -1

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-1,-0.5,0,0.5,1), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-1,ub=1)

## c = -1.5 -----------------------------------------------------------------
cT <- -1.5

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-1,-0.5,0,0.5,1), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-1,ub=1)

## c = -2 -------------------------------------------------------------------
cT <- -2

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-1,-0.5,0,0.5,1), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-1,ub=1)


#############################################################################
## Distribution 4
## Test density with pole.
#############################################################################

lf <- function(x) { -log(abs(x))/2 }    ## 1/sqrt(x)
dlf <- function(x) { -1/(2*x) }
## extrema: 0
d2lf <- function(x) { 1/(2*x^2) }

## c = -1.5 -----------------------------------------------------------------
cT <- -1.5

gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-1,0,1), cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-1,ub=1)


#############################################################################
## Distribution 5
## Test different values for 'c'
#############################################################################

lf <- function(x) { -2*x^4 + 4*x^2 } 
dlf <- function(x) { -8*x^3 + 8*x }
d2lf <- function(x) { -24*x^2+8 }
## extrema: -1, 0, 1

cT <- c(-0.5, 2, -2, 0.5, -1, 0)
ib <- c(-Inf,-2, -1,   0,  1, 2, Inf)

gen <- Tinflex.setup(lf, dlf, d2lf, ib=ib, cT=cT, rho=rho)
run.test(gen,cT,lf)

cT <- c(-0.5, 2, -2, 0.5, -1, 0)
ib <- c(  -3,-2, -1,   0,  1, 2, 3)

gen <- Tinflex.setup(lf, dlf, d2lf, ib=ib, cT=cT, rho=rho)
run.test(gen,cT,lf,lb=-3,ub=3)

rm(cT, ib)


#############################################################################
## Distribution 6
## Test density with infection point at interval boundary
#############################################################################

## c = 0 --------------------------------------------------------------------
cT <- 0

lf <- function(x) { -x^4+6*x^2 }
dlf <- function(x) { 12*x-4*x^3 }
d2lf <- function(x) { 12-12*x^2 }

## inflection points: -1, 1
gen <- Tinflex.setup(lf, dlf, d2lf, ib=c(-Inf,-2,-1,0,1,2,Inf), cT=cT, rho=rho)
print(gen)
run.test(gen,cT,lf)


#############################################################################
## Summary of tests.
#############################################################################

## Number of tests.
length(pvals)

## p-values for goodness-of-fit tests.
summary(pvals)
  
## Level-2 test for p-values (must be uniformly distributed).
if (length(pvals)>0)
  ks.test(x=pvals, y="punif", alternative = "greater")

## Print number of tests where R and C code return (slightly) different results.
comp.fail

## Stop timer.
run.time <- (proc.time() - time.start)[3]  ## "elapsed" time
run.time

## Failed tests?
if (length(pvals[pvals < pval.threshold])) {
  stop (paste(length(pvals[pvals < pval.threshold]),"out of",length(pvals),
              "goodness-of-fit tests failed!"))
} else {
  cat ("All goodness-of-fit tests passed!\n")
}
  

#############################################################################
