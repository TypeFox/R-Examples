#############################################################################
##
##  List of Generation Methods
##
#############################################################################

## --------------------------------------------------------------------------
## Type of method
## Important: The spelling of these words is crucial

METHODS.TYPE <- c("Automatic    ","Inversion    ","Rejection    ","Select method")

## --------------------------------------------------------------------------
## Methods for continuous univariate distributions

METHODS.CONT <-
  array(c(##  Name                                             ## Symbol
          "-- select a method for continous distributions --", NA,
          "ARS -- Adaptive Rejection Sampling",                "ARS",
          "ITDR -- Inverse Transformed Density Rejection",     "ITDR",
          "PINV -- Polynomial interpolation of INVerse CDF",   "PINV",
          "SROU -- Simple Ratio-Of-Uniforms Method",           "SROU",
          "TDR -- Transformed Density Rejection",              "TDR"
          ))

METHODS.CONT.AUTO <-
  array(c(##  Name                                             ## Symbol
          "Automatic method",                                  "PINV",                                 
          "Inversion method",                                  "PINV",                                 
          "Rejection method",                                  "TDR"                                 
          ))

## --------------------------------------------------------------------------
## Methods for discrete univariate distributions

METHODS.DISCR <-
  array(c(##  Name                                             ## Symbol
          "-- select a method for discrete distributions --",  NA,
          "DAU -- Alias-Urn Method",                           "DAU",
          "DGT -- Table guided discrete inversion",            "DGT",
          "DARI -- Discrete Automatic Rejection Inversion",    "DARI"
          ))

METHODS.DISCR.AUTO <-
  array(c(##  Name                                             ## Symbol
          "Automatic method",                                  "DGT",                                 
          "Inversion method",                                  "DGT",                                 
          "Rejection method",                                  "DARI"                                 
          ))

## --------------------------------------------------------------------------

## combine in single list
METHODS <- list(continuous=METHODS.CONT, continuous.auto=METHODS.CONT.AUTO,
                discrete=METHODS.DISCR,  discrete.auto  =METHODS.DISCR.AUTO)

## set dimensions and rownames
for (l in names(METHODS)) {
  dim(METHODS[[l]]) <- c(2, length(METHODS[[l]])/2)
  rownames(METHODS[[l]]) <- c("name","symbol")
}


#############################################################################
##
##  Widgets for inserting parameters for generation methods
##
#############################################################################

## --------------------------------------------------------------------------
## Synopsis:
##
##   param.distr.METH ( main, group, args )
##
##   METH ... symbol for method (see above lists) in upper case letters 
##
##   Arguments:
##     main  ... main window
##     group ... container where widgets for inserting parameters are placed
##     args  ... a _list_ of arguments for the consructor of the generator
##               object. The 'names' attribute of 'args' gives the argument
##               its values are the defaults for the corresponding argument.
##
## Important!
##  There _must_ such a function for _each_ generation method in the above lists.
## --------------------------------------------------------------------------

## --------------------------------------------------------------------------
## generic function

param.method.generic <- function(main, group) {
  ## frame for parameters
  frame <- gframe(" Parameters ",
                  horizontal=FALSE, container=group)
  ## there are no parameters
  glabel(" -- none -- ", container=frame, anchor=c(-1,0))

  ## return empty list
  return(list())
}


## --------------------------------------------------------------------------
## Methods for continuous univariate distributions
## --------------------------------------------------------------------------

## --------------------------------------------------------------------------
## ARS

param.method.ARS <- function(main, group) {
  param.method.generic(main, group)
}


## --------------------------------------------------------------------------
## ITDR

param.method.ITDR <- function(main, group) {
  param.method.generic(main, group)
}


## --------------------------------------------------------------------------
## PINV

param.method.PINV <- function(main,group) {

  ## frame for parameters
  frame <- gframe(" Parameters ", container=group)

  ## parse arguments of constructor
  args <- formals( tag(main,"method")[["constr"]] )

  ## store widgets for parameters in a list
  params.gwl <- list()

  ## use a tabular layout for widgets
  tbl <- glayout(container=frame)

  params.gwl[["uresolution"]] <-
    gedit(text=args[["uresolution"]], coerce.with=as.numeric, width=15, container=tbl)
  tooltip(params.gwl[["uresolution"]]) <-
    "Maximal tolerated error in uniform scale (u-error)"
  tag(params.gwl[["uresolution"]],"label") <- "u-resolution"
  tbl[1,2,anchor=c(-1,0)] <- "u-resolution"
  tbl[1,3,anchor=c(-1,0)] <- params.gwl[["uresolution"]]

  params.gwl[["smooth"]] <-
    gcheckbox("smooth inverse function", checked=FALSE, container=tbl)
  tooltip(params.gwl[["smooth"]]) <-
    "Check this box if approximate inverse CDF must be differentiable.\n(Not recommended)"
  tbl[2,2:3,anchor=c(-1,0)] <- params.gwl[["smooth"]]

  return(params.gwl)
}


## --------------------------------------------------------------------------
## SROU

param.method.SROU <- function(main,group) {

  ## frame for parameters
  frame <- gframe(" Parameters ", container=group)

  ## parse arguments of constructor
  args <- formals( tag(main,"method")[["constr"]] )

  ## store widgets for parameters in a list
  params.gwl <- list()

  ## use a tabular layout for widgets
  tbl <- glayout(container=frame)

  params.gwl[["r"]] <-
    gedit(text=args[["r"]],coerce.with=as.numeric,container=tbl)
  tooltip(params.gwl[["r"]]) <-
    paste("Use higher values for distributions with heavy tails",
          "(see help page for details).", sep="\n")

  tbl[1,2,anchor=c(-1,0)] <- "concavity parameter r"
  tbl[1,3,anchor=c(-1,0)] <- params.gwl[["r"]]

  return(params.gwl)
}


## --------------------------------------------------------------------------
## TDR

param.method.TDR <- function(main, group) {
  param.method.generic(main, group)
}


## --------------------------------------------------------------------------
## Methods for discrete univariate distributions
## --------------------------------------------------------------------------

## --------------------------------------------------------------------------
## DARI

param.method.DARI <- function(main, group) {
  param.method.generic(main, group)
}


## --------------------------------------------------------------------------
## DAU

param.method.DAU <- function(main, group) {
  param.method.generic(main, group)
}


## --------------------------------------------------------------------------
## DGT

param.method.DGT <- function(main, group) {
  param.method.generic(main, group)
}

## --------------------------------------------------------------------------


#############################################################################
##
##  Widgets for inserting parameters for user-defined methods
##
#############################################################################

## --------------------------------------------------------------------------
## Methods for continuous univariate distributions
## --------------------------------------------------------------------------

## --------------------------------------------------------------------------
## ARS

param.distr.ARS <- function(main, group) {
  ## synopsis: ars.new(logpdf, dlogpdf=NULL, lb, ub, ...)

  ## store widgets for parameters in a list
  params.gwl <- list()

  ## use a tabular layout for widgets
  tbl <- glayout(container=group)

  ## density function (require --> return NA in case of an empty field)
  params.gwl[["pdf"]] <- gedit(text="", coerce.with=wrap.function.body.NA, container=tbl)
  tooltip(params.gwl[["pdf"]]) <-
    paste("Log-density of distribution.",
          "Do not forget to set domain of distribution below!", sep="\n")
  tag(params.gwl[["pdf"]],"label") <- "log-density"
  tbl[1,2,anchor=c(-1,0)] <- "log-density"
  tbl[1,3,anchor=c(-1,0)] <- "<- function(x) {"
  tbl[1,4,anchor=c(-1,0)] <- params.gwl[["pdf"]]
  tbl[1,5,anchor=c(-1,0)] <- "}"

  params.gwl[["dpdf"]] <- gedit(text="", coerce.with=wrap.function.body.NA, container=tbl)
  tooltip(params.gwl[["dpdf"]]) <- "Derivative of given log-density"
  tag(params.gwl[["dpdf"]],"label") <- "derivative"
  tbl[2,2,anchor=c(-1,0)] <- "derivative"
  tbl[2,3,anchor=c(-1,0)] <- "<- function(x) {"
  tbl[2,4,anchor=c(-1,0)] <- params.gwl[["dpdf"]]
  tbl[2,5,anchor=c(-1,0)] <- "}"
  
  ## we need 'islog=TRUE': show disabled checkbox 
  params.gwl[["islog"]] <- 
    gcheckbox("log-density", checked=TRUE, container=tbl)
  enabled(params.gwl[["islog"]]) <- FALSE
  tbl[3,2:3,anchor=c(-1,0)] <- params.gwl[["islog"]]
  
  return(params.gwl)
}


## --------------------------------------------------------------------------
## ITDR

param.distr.ITDR <- function(main, group) {
  ## synopsis: itdr.new(pdf, dpdf, lb, ub, pole, islog=FALSE, ...)

  ## store widgets for parameters in a list
  params.gwl <- list()
  
  ## use a tabular layout for widgets
  tbl <- glayout(container=group)

  ## density function (require --> return NA in case of an empty field)
  params.gwl[["pdf"]] <- gedit(text="", coerce.with=wrap.function.body.NA, container=tbl)
  tooltip(params.gwl[["pdf"]]) <-
    paste("(Log-) density of distribution.",
          "Do not forget to set domain of distribution below!", sep="\n")
  tag(params.gwl[["pdf"]],"label") <- "density"
  tbl[1,2,  anchor=c(-1,0)] <- "density"
  tbl[1,3,  anchor=c(-1,0)] <- "<- function(x) {"
  tbl[1,4:5,anchor=c(-1,0)] <- params.gwl[["pdf"]]
  tbl[1,6,  anchor=c(-1,0)] <- "}"

  params.gwl[["dpdf"]] <- gedit(text="", coerce.with=wrap.function.body.NA, container=tbl)
  tooltip(params.gwl[["dpdf"]]) <- "Derivative of given (log-) density"
  tag(params.gwl[["dpdf"]],"label") <- "derivative"
  tbl[2,2,  anchor=c(-1,0)] <- "derivative"
  tbl[2,3,  anchor=c(-1,0)] <- "<- function(x) {"
  tbl[2,4:5,anchor=c(-1,0)] <- params.gwl[["dpdf"]]
  tbl[2,6,  anchor=c(-1,0)] <- "}"
  
  ## do we have log-density?
  params.gwl[["islog"]] <- 
    gcheckbox("log-density", checked=FALSE, container=tbl)
  tooltip(params.gwl[["islog"]]) <-
    "Check this box if the given expressions are the log-density and its derivative."
  tbl[3,2:3,anchor=c(-1,0)] <- params.gwl[["islog"]]
  
  ## pole (it is stored in the distribution object using argument 'mode')
  params.gwl[["mode"]] <-
    gedit(text="", coerce.with=as.numeric, width=15, container=tbl)
  tooltip(params.gwl[["mode"]]) <-
    "Set pole of density function"
  tag(params.gwl[["mode"]],"label") <- "pole of distribution"
  tbl[4,2:3,anchor=c(-1,0)] <- "pole of distribution ('mode')"
  tbl[4,4,  anchor=c(-1,0)] <- params.gwl[["mode"]]
  
  return(params.gwl)
}


## --------------------------------------------------------------------------
## PINV

param.distr.PINV <- function(main, group) {
  ## synopsis: pinv.new(pdf, cdf, lb, ub, islog=FALSE, center=0,
  ##                    uresolution=1.e-10, smooth=FALSE, ...)

  ## store widgets for parameters in a list
  params.gwl <- list()

  ## use a tabular layout for widgets
  tbl <- glayout(container=group)

  tbl[1,2:3,anchor=c(-1,0)] <- "Need one of ..."
  
  ## density function (my be empty --> return NULL in case of an empty field)
  params.gwl[["pdf"]] <- gedit(text="", coerce.with=wrap.function.body.NULL, width=50, container=tbl)
  tooltip(params.gwl[["pdf"]]) <-
    paste("(Log-) density of distribution.",
          "Do not forget to set domain of distribution below!", sep="\n")
  tag(params.gwl[["pdf"]],"label") <- "density"
  tbl[2,2,  anchor=c(-1,0)] <- "density"
  tbl[2,3,  anchor=c(-1,0)] <- "<- function(x) {"
  tbl[2,4:5,anchor=c(-1,0)] <- params.gwl[["pdf"]]
  tbl[2,6,  anchor=c(-1,0)] <- "}"

  ## distribution function (my be empty --> return NULL in case of an empty field)
  params.gwl[["cdf"]] <- gedit(text="", coerce.with=wrap.function.body.NULL, width=50, container=tbl)
  tooltip(params.gwl[["cdf"]]) <-
    paste("(Logarithm of) distribution function.",
          "Only used if density is not given!", sep="\n")
  ## tag(params.gwl[["cdf"]],"label") <- "cdf"
  tbl[3,2,  anchor=c(-1,0)] <- "cdf"
  tbl[3,3,  anchor=c(-1,0)] <- "<- function(x) {"
  tbl[3,4:5,anchor=c(-1,0)] <- params.gwl[["cdf"]]
  tbl[3,6,  anchor=c(-1,0)] <- "}"

  ## do we have log-density?
  params.gwl[["islog"]] <- 
    gcheckbox("log-density", checked=FALSE, container=tbl)
  tooltip(params.gwl[["islog"]]) <-
    "Check this box if the given expressions are the logarithms of the respective functions."
  tbl[4,2:3,anchor=c(-1,0)] <- params.gwl[["islog"]]
  
  ## center of distribution
  params.gwl[["center"]] <- 
    gedit(text="0", coerce.with=as.numeric, width=15, container=tbl)
  tooltip(params.gwl[["center"]]) <-
    paste("Center of distribution",
          "A 'typical' point of the distribution.",
          "Give a point near the mode the distribution.",
          sep="\n")
  tbl[5,2:3,anchor=c(-1,0)] <- "center of distribution"
  tbl[5,4,  anchor=c(-1,0)] <- params.gwl[["center"]]

  return(params.gwl)
}


## --------------------------------------------------------------------------
## SROU

param.distr.SROU <- function(main, group) {
  ## synopsis: srou.new(pdf, lb, ub, mode, area, islog=FALSE, r=1, ...)

  ## store widgets for parameters in a list
  params.gwl <- list()

  ## use a tabular layout for widgets
  tbl <- glayout(container=group)

  ## density function (require --> return NA in case of an empty field)
  params.gwl[["pdf"]] <- gedit(text="", coerce.with=wrap.function.body.NA, container=tbl)
  tooltip(params.gwl[["pdf"]]) <-
    paste("(Log-) density of distribution.",
          "Do not forget to set domain of distribution below!", sep="\n")
  tag(params.gwl[["pdf"]],"label") <- "density"
  tbl[1,2,  anchor=c(-1,0)] <- "density"
  tbl[1,3,  anchor=c(-1,0)] <- "<- function(x) {"
  tbl[1,4:5,anchor=c(-1,0)] <- params.gwl[["pdf"]]
  tbl[1,6,  anchor=c(-1,0)] <- "}"

  ## do we have log-density?
  params.gwl[["islog"]] <- 
    gcheckbox("log-density", checked=FALSE, container=tbl)
  tooltip(params.gwl[["islog"]]) <-
    "Check this box if the given expression is the log-density."
  tbl[2,2:3,anchor=c(-1,0)] <- params.gwl[["islog"]]
  
  ## mode
  params.gwl[["mode"]] <- 
    gedit(text="", coerce.with=as.numeric, width=15, container=tbl)
  tooltip(params.gwl[["mode"]]) <- "Mode of distribution"
  tag(params.gwl[["mode"]],"label") <- "mode of distribution"
  tbl[3,2:3,anchor=c(-1,0)] <- "mode of distribution"
  tbl[3,4,anchor=c(-1,0)] <- params.gwl[["mode"]]

  ## area below density
  params.gwl[["area"]] <- 
    gedit(text="", coerce.with=as.numeric, width=15, container=tbl)
  tooltip(params.gwl[["area"]]) <- "Area below density"
  tag(params.gwl[["area"]],"label") <- "area below density"
  tbl[4,2:3,anchor=c(-1,0)] <- "area below density"
  tbl[4,4,anchor=c(-1,0)] <- params.gwl[["area"]]
  
  return(params.gwl)
}


## --------------------------------------------------------------------------
## TDR

param.distr.TDR <- function(main, group) {
  ## synopsis: tdr.new(pdf, dpdf=NULL, lb, ub, islog=FALSE, ...)

  ## store widgets for parameters in a list
  params.gwl <- list()

  ## use a tabular layout for widgets
  tbl <- glayout(container=group)

  ## density function (require --> return NA in case of an empty field)
  params.gwl[["pdf"]] <- gedit(text="", coerce.with=wrap.function.body.NA, container=tbl)
  tooltip(params.gwl[["pdf"]]) <-
    paste("(Log-) density of distribution.",
          "Do not forget to set domain of distribution below!", sep="\n")
  tag(params.gwl[["pdf"]],"label") <- "density"
  tbl[1,2,anchor=c(-1,0)] <- "density"
  tbl[1,3,anchor=c(-1,0)] <- "<- function(x) {"
  tbl[1,4,anchor=c(-1,0)] <- params.gwl[["pdf"]]
  tbl[1,5,anchor=c(-1,0)] <- "}"

  params.gwl[["dpdf"]] <- gedit(text="", coerce.with=wrap.function.body.NA, container=tbl)
  tooltip(params.gwl[["dpdf"]]) <- "Derivative of given (log-) density"
  tag(params.gwl[["dpdf"]],"label") <- "derivative"
  tbl[2,2,anchor=c(-1,0)] <- "derivative"
  tbl[2,3,anchor=c(-1,0)] <- "<- function(x) {"
  tbl[2,4,anchor=c(-1,0)] <- params.gwl[["dpdf"]]
  tbl[2,5,anchor=c(-1,0)] <- "}"
  
  ## do we have log-density?
  params.gwl[["islog"]] <- 
    gcheckbox("log-density", checked=FALSE, container=tbl)
  tooltip(params.gwl[["islog"]]) <-
    "Check this box if the given expressions are the log-density and its derivative."
  tbl[3,2:3,anchor=c(-1,0)] <- params.gwl[["islog"]]
  
  return(params.gwl)
}


## --------------------------------------------------------------------------
## Methods for discrete univariate distributions
## --------------------------------------------------------------------------

## --------------------------------------------------------------------------
## DARI

param.distr.DARI <- function(main, group) {
  ## synopsis: dari.new(pmf, lb, ub, mode=NA, sum=1, ...)

  ## store widgets for parameters in a list
  params.gwl <- list()

  ## use a tabular layout for widgets
  tbl <- glayout(container=group)

  ## probability function (require --> return NA in case of an empty field)
  params.gwl[["pmf"]] <- gedit(text="", coerce.with=wrap.function.body.NA, container=tbl)
  tooltip(params.gwl[["pmf"]]) <-
    paste("probability mass function of distribution.",
          "Do not forget to set domain of distribution below!", sep="\n")
  tag(params.gwl[["pmf"]],"label") <- "pmf"
  tbl[1,2,anchor=c(-1,0)] <- "pmf"
  tbl[1,3,anchor=c(-1,0)] <- "<- function(x) {"
  tbl[1,4,anchor=c(-1,0)] <- params.gwl[["pmf"]]
  tbl[1,5,anchor=c(-1,0)] <- "}"

  ## mode
  params.gwl[["mode"]] <- 
    gedit(text="", coerce.with=as.integer, width=15, container=tbl)
  tooltip(params.gwl[["mode"]]) <- "Mode of distribution"
  tag(params.gwl[["mode"]],"label") <- "mode of distribution"
  tbl[3,2:3,anchor=c(-1,0)] <- "mode of distribution"
  tbl[3,4,anchor=c(-1,0)] <- params.gwl[["mode"]]

  ## sum of probabilities
  params.gwl[["sum"]] <- 
    gedit(text="", coerce.with=as.numeric, width=15, container=tbl)
  tooltip(params.gwl[["sum"]]) <- "Sum of probabilities"
  tag(params.gwl[["sum"]],"label") <- "sum of probablities"
  tbl[4,2:3,anchor=c(-1,0)] <- "sum of probabilities"
  tbl[4,4,anchor=c(-1,0)] <- params.gwl[["sum"]]
  tbl[4,5,anchor=c(-1,0)] <- "(approx.)"
  
  return(params.gwl)
}


## --------------------------------------------------------------------------
## DAU

param.distr.DAU <- function(main, group) {
  ## synopsis: dau.new(pv, from=1)

  ## store widgets for parameters in a list
  params.gwl <- list()

  ## use a tabular layout for widgets
  tbl <- glayout(container=group)

  ## density function (require --> return NA in case of an empty field)
  params.gwl[["pv"]] <- gedit(text="", container=tbl)
  tooltip(params.gwl[["pv"]]) <-
    paste("A vector of probability weights.",
          "(Need not sum to one, but non-negative and not all zero.)", sep="\n")
  tag(params.gwl[["pv"]],"label") <- "probability vector"
  tbl[1,2,anchor=c(-1,0)] <- "probability vector"
  tbl[1,3,anchor=c(-1,0)] <- "<-"
  tbl[1,4,anchor=c(-1,0)] <- params.gwl[["pv"]]
  
  return(params.gwl)
}


## --------------------------------------------------------------------------
## DGT

param.distr.DGT <- function(main, group) {
  ## synopsis: dgt.new(pv, from=1)
  ## Same as dau.new.
  param.distr.DAU(main, group)
}

## --------------------------------------------------------------------------
