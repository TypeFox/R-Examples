#source("d:/MWP/eAnalysis/regtest/R/regtest.R")

#! \name{timefactor}
#! \alias{timefactor}
#! \title{ compare timing of two expressions }
#! \description{
#!   Compare timings of two expressions, expressions that do not take enough time to be measured can be repeated often enough to be measured.
#! }
#! \usage{
#! timefactor(nom, denom, repnom = 1, repdenom = 1)
#! }
#! \arguments{
#!   \item{nom}{ nominator expression }
#!   \item{denom}{ denominator expression }
#!   \item{repnom}{ no. of repetitions of nominator }
#!   \item{repdenom}{ no. of repetitions of denominator }
#! }
#! \value{
#!   matrix with absolzute and relative timings
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{system.time}}, \code{\link{Sys.sleep}} }
#! \examples{
#!   timefactor(Sys.sleep(0.1), Sys.sleep(1), 10, 1)
#! }
#! \keyword{ misc }


timefactor <-
function (nom, denom, repnom = 1, repdenom = 1)
{
    gc()
    n <- system.time(for (f in 1:repnom) eval(substitute(nom), envir=parent.frame()))
    d <- system.time(for (f in 1:repdenom) eval(substitute(denom), envir=parent.frame()))
    cbind(nom = n/repnom, denom = d/repdenom, factor = repdenom/repnom *
        n/d)[c(1:3), ]
}




#! \name{is.all.equal}
#! \alias{is.all.equal}
#! \title{ wrapper for all.equal }
#! \description{
#!   like all.equal but always returns logical
#! }
#! \usage{
#! is.all.equal(a, b, \dots)
#! }
#! \arguments{
#!   \item{a}{ expression to compare }
#!   \item{b}{ expression to compare }
#!   \item{\dots}{ further arguments to \code{\link{all.equal}} }
#! }
#! \value{ TRUE or FALSE }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{all.equal}}, \code{\link{identical}}, \code{\link{binregtest}} }
#! \examples{
#!   all.equal(1,2)
#!   is.all.equal(1,2)
#! }
#! \keyword{ debugging }
#! \keyword{ utilities }

is.all.equal <- function(a, b, ...){
  x <- all.equal(a,b,...)
  is.logical(x) && x
}


#! \name{binregtest}
#! \alias{binregtest}
#! \title{ Binary regression test }
#! \description{
#!   This function compares to parallel implementations of the same function for equality with respect to all (defined) parameter combinations.
#! }
#! \usage{
#! binregtest(FUN1, FUN2, ..., PARS = NULL, PAR1 = NULL, PAR2 = NULL, WHICH = sample(1:n), TRYALL = TRUE, COMP = is.all.equal, NAME = "UNNAMED binregtest", VERBOSE = FALSE)
#! }
#! \arguments{
#!   \item{FUN1}{ first function }
#!   \item{FUN2}{ second function }
#!   \item{\dots}{ common arguments, each specified as a list }
#!   \item{PARS}{ list of common arguments, each specified as a list, helps to use argument names that are used by \code{binregtest} itself}
#!   \item{PAR1}{ optional parameters only handed over to FUN1 (default NULL) }
#!   \item{PAR2}{ optional parameters only handed over to FUN2 (default NULL) }
#!   \item{WHICH}{ optional integer vector defining a subset of the possible parameter combinations }
#!   \item{TRYALL}{ FALSE to interrupt testing once an error has been found (default TRUE tests everything) }
#!   \item{COMP}{ function to compare results (default \code{\link{is.all.equal}}) }
#!   \item{NAME}{ character scalar describing this regression test }
#!   \item{VERBOSE}{ TRUE to verbose all tests (default FALSE) }
#! }
#! \value{
#!   TRUE if all tests were successful, FALSE otherwise
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{try}}, \code{\link{is.all.equal}} }
#! \examples{
#! wronglog <- function(x, base=exp(1)){
#!   if (x>0)
#!     log(x, base=base)
#!   else
#!     NA
#! }
#! binregtest(wronglog, log, x=as.list(0:3), base=list(2, exp(1), 10))
#! }
#! \keyword{ debugging }
#! \keyword{ documentation }


binregtest <- function(
  FUN1
, FUN2
, ...
, PARS  = NULL
, PAR1  = NULL
, PAR2  = NULL
, WHICH = sample(1:n)
, TRYALL= TRUE
, COMP  = is.all.equal
, NAME  = "UNNAMED binregtest"
, VERBOSE = FALSE
){
  FUN1NAME <- deparse(substitute(FUN1))
  FUN2NAME <- deparse(substitute(FUN2))
  PAR1NAME <- deparse(substitute(PAR1))
  PAR2NAME <- deparse(substitute(PAR2))
  all(sapply(list(...), is.list))
  all(sapply(PARS, is.list))
  PARS <- c(list(...), PARS)
  np <- length(PARS)
  ip <- as.integer(seq(np))
  names(ip) <- names(PARS)
  N <- sapply(PARS, length)
  n <- rev(cumprod(N))[1]
  L <- lapply(PARS, function(i)as.integer(seq(along=i)))
  success <- TRUE
  d <- do.call("expand.grid", c(L, list(KEEP.OUT.ATTRS = FALSE)))
  for (i in WHICH){
    di <- d[i,,drop=TRUE]
    P <- lapply(ip, function(j){
      PARS[[j]][[di[[j]]]]
    })
    if (VERBOSE){
      cat("\n\nVERBOSEing ", NAME," ", i, " of ", length(WHICH),"\n", sep="")
      str(FUN1)
      str(PAR1)
      str(FUN2)
      str(PAR2)
      str(P)
    }
      str
    for (j in rev(ip))
      if (identical(P[[j]], missing))
        P[[j]] <- NULL
    seed <- as.integer(runif(1, min=0, max=.Machine$integer.max))
    set.seed(seed)
    RET1 <- try(do.call("FUN1", c(PAR1, P)))
    set.seed(seed)
    RET2 <- try(do.call("FUN2", c(PAR2, P)))
    this.success <- TRUE
    if (inherits(RET1, "try-error")){
      this.success <- FALSE
      cat("\n\n",NAME,"\n",sep="")
      cat("error in FUN1", FUN1NAME, "PAR1", PAR1NAME, "\n")
      str(FUN1)
      str(PAR1)
    }
    if (inherits(RET2, "try-error")){
      this.success <- FALSE
      cat("\n\n",NAME,"\n",sep="")
      cat("error in FUN2", FUN2NAME, "PAR2", PAR2NAME, "\n")
      str(FUN2)
      str(PAR2)
    }
    if (!this.success){
      cat("common PARS\n")
      str(P)
    }else{
      if (!COMP(RET1, RET2)){
        this.success <- FALSE
        cat("\n\n",NAME,"\n",sep="")
        cat("difference between FUN1", FUN1NAME, "PAR1", PAR1NAME, "and FUN2", FUN2NAME, "PAR2", PAR2NAME, "\n")
        str(FUN1)
        str(PAR1)
        str(FUN2)
        str(PAR2)
        cat("common PARS\n")
        str(P)
        cat("RET1\n")
        str(RET1)
        cat("RET2\n")
        str(RET2)
      }
    }
    success <- success && this.success
    if (!this.success && !TRYALL)
      break
  }
  success
}
