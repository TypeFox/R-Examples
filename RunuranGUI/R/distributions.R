#############################################################################
##
##  List of Built-in Distributions
##
#############################################################################

## --------------------------------------------------------------------------
## Type of distribution
## Important: The spelling of these words is crucial

DISTRIBUTIONS.TYPE <- c("continuous","discrete")

## --------------------------------------------------------------------------
## How distribution is defined
## Important: The spelling of these words is crucial

DISTRIBUTIONS.HOWDEF <- c("built-in","user-defined")

## --------------------------------------------------------------------------
## Continuous univariate distributions

DISTRIBUTIONS.CONT <-
  array(c(##  Name                                    ## Symbol
          "-- select a continuous distribution --",   NA,
          "Beta distribution",                        "beta",
          "Cauchy distribution",                      "cauchy",
          "Chi distribution",                         "chi",      ### qchi
          "Chi-squared distribution",                 "chisq",
          "Exponential distribution",                 "exp",
          "F distribution",                           "f",
          "Frechet (extreme value II) distribution",  "frechet",
          "Gamma distribution",                       "gamma",
          "Generalized hyperbolic distribution",      "ghyp",     ### qghyp
          "Generalized inverse Gaussian distribution","gig",      ### qgig
          "Gumbel (extreme value I) distribution",    "gumbel",   ### qgumbel
          "Hyperbolic distribution",                  "hyperbolic", # ghyperbolic
          "Inverse Gaussian distribution",            "ig",       ### qig
          "Laplace (double exponential) distribution","laplace",  ### qlaplace
          "Lognormal distribution",                   "lnorm",
          "Logistic distribution",                    "logis",
          "Lomax (Pareto II) distribution",           "lomax",    ### qlomax
          "Normal distribution",                      "norm",
          "Powerexponential (Subbotin) distribution", "powerexp", ### qpowerexp
          "Rayleigh distribution",                    "rayleigh", ### qrayleigh
          "Slash distribution",                       "slash",
          "Student's t distribution",                 "t",
          "Weibull distribution",                     "weibull"
          ))

## --------------------------------------------------------------------------
## Discrete univariate distributions

DISTRIBUTIONS.DISCR <-
  array(c(##  Name                                ## Symbol
          "-- select a discrete distribution --", NA,
          "Binomial distribution",                "binom",
          "Geometric distribution",               "geom",
          "Hypergeometric distribution",          "hyper",
          "Logarithmic distribution",             "logarithmic",
          "Negative binomial distribution",       "nbinom",
          "Poisson distribution",                 "pois"
          ))


## --------------------------------------------------------------------------

## combine in single list
DISTRIBUTIONS <- list(continuous=DISTRIBUTIONS.CONT, discrete=DISTRIBUTIONS.DISCR)

## set dimensions and rownames
for (l in names(DISTRIBUTIONS)) {
  dim(DISTRIBUTIONS[[l]]) <- c(2, length(DISTRIBUTIONS[[l]])/2)
  rownames(DISTRIBUTIONS[[l]]) <- c("name","symbol")
}


#############################################################################
##
##  Widgets for inserting parameters for built-in distributions
##
#############################################################################

## --------------------------------------------------------------------------
## Synopsis:
##
##   param.distr.DISTR ( main, group )
##
##   DISTR ... symbol for distribution (see above lists) in upper case letters 
##
##   Arguments:
##     main  ... main window
##     group ... container where widgets for inserting parameters are placed
##
## Remark:
##  There is no need to provide such a function for each distribution in the
##  above lists. If such a function is not presented the widget is created
##  using argument names and defaults from the corresponding function
##  "ud...".
## --------------------------------------------------------------------------

## --------------------------------------------------------------------------
## generic function

param.distr.generic <- function(main,group) {

  ## get data about distribution
  constr <- tag(main,"distribution")$constr

  ## store widgets for parameters in a list
  params.gwl <- list()

  ## parse arguments of constructor
  args <- function.args(constr)

  ## remove arguments for domain
  args <- args[-match(c("lb","ub"), names(args))]

  ## create widgets (use a tabular layout)
  tbl <- glayout(container=group)

  if (length(args)>0) {
    for (i in 1:length(args)) {
      p <- names(args)[i]
      params.gwl[[p]] <- gedit(text=args[[p]], width=20,
                               coerce.with=as.numeric, container=tbl)
      tbl[i,1,anchor=c(1,0)] <- p
      tbl[i,2] <- params.gwl[[p]]
    }
  }
  else {
    tbl[1,1,anchor=c(1,0)] <- "-- none --"
    tbl[1,2,anchor=c(1,0)] <- ""
  }    

  return(params.gwl)
}


## --------------------------------------------------------------------------
## parameters for distribution NORM
##
##EXPERIMENTAL.param.distr.NORM <- function(main,group) {
##
##  ## get data about distribution
##  constr <- tag(main,"distribution")$constr
##
##  ## store widgets for parameters in a list
##  params.gwl <- list()
##
##  ## parse arguments of constructor
##  args <- function.args(constr)
##
##  ## create widgets (use a tabular layout)
##  tbl <- glayout(container=group)
##
##  params.gwl[["mean"]] <- gedit(text=args[["mean"]], coerce.with=as.numeric, container=tbl)
##  tbl[1,1,anchor=c(1,0)] <- "mean"
##  tbl[1,2] <- params.gwl[["mean"]]
##
##  params.gwl[["sd"]] <- gedit(text=args[["sd"]], coerce.with=as.numeric, container=tbl)
##  tbl[2,1,anchor=c(1,0)] <- "sd"
##  tbl[2,2] <- params.gwl[["sd"]]
##  tbl[2,3,anchor=c(-1,0)] <- "standard deviation (positive)"
##
####    addHandlerBlur(vars[[p]], function(h, ...) {
####      str <- svalue(h$obj)
####      if (nchar(str)==4)
####        font(h$obj) <- c(color="red")
####      else
####        font(h$obj) <- c(color="black")
####    })
##
##  return(params.gwl)
##}
##
## --------------------------------------------------------------------------

