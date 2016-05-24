#' Simulate EU-SILC population data
#' 
#' Simulate population data for the European Statistics on Income and Living
#' Conditions (EU-SILC).
#' 
#' @name simEUSILC
#' @param dataS a \code{data.frame} containing EU-SILC survey data.
#' @param hid a character string specifying the column of \code{dataS} that
#' contains the household ID.
#' @param wh a character string specifying the column of \code{dataS} that
#' contains the household sample weights.
#' @param wp a character string specifying the column of \code{dataS} that
#' contains the personal sample weights.
#' @param hsize an optional character string specifying a column of
#' \code{dataS} that contains the household size. If \code{NULL}, the household
#' sizes are computed.
#' @param strata a character string specifying the column of \code{dataS} that
#' define strata. Note that this is currently a required argument and only one
#' stratification variable is supported.
#' @param pid an optional character string specifying a column of \code{dataS}
#' that contains the personal ID.
#' @param age a character string specifying the column of \code{dataS} that
#' contains the age of the persons (to be used for setting up the household
#' structure).
#' @param gender a character string specifying the column of \code{dataS} that
#' contains the gender of the persons (to be used for setting up the household
#' structure).
#' @param categorizeAge a logical indicating whether age categories should be
#' used for simulating additional categorical and continuous variables to
#' decrease computation time.
#' @param breaksAge numeric; if \code{categorizeAge} is \code{TRUE}, an
#' optional vector of two or more break points for constructing age categories,
#' otherwise ignored.
#' @param categorical a character vector specifying additional categorical
#' variables of \code{dataS} that should be simulated for the population data.
#' @param income a character string specifying the variable of \code{dataS}
#' that contains the personal income (to be simulated for the population data).
#' @param method a character string specifying the method to be used for
#' simulating personal income. Accepted values are \code{"multinom"} (for using
#' multinomial log-linear models combined with random draws from the resulting
#' ategories) and \code{"twostep"} (for using two-step regression models
#' combined with random error terms).
#' @param breaks if \code{method} is \code{"multinom"}, an optional numeric
#' vector of two or more break points for categorizing the personal income. If
#' missing, break points are computed using weighted quantiles.
#' @param lower,upper numeric values; if \code{method} is \code{"multinom"} and
#' \code{breaks} is \code{NULL}, these can be used to specify lower and upper
#' bounds other than minimum and maximum, respectively. Note that if \code{gpd}
#' is \code{TRUE} (see below), \code{upper} defaults to \code{Inf}.
#' @param equidist logical; if \code{method} is \code{"multinom"} and
#' \code{breaks} is \code{NULL}, this indicates whether the (positive) default
#' break points should be equidistant or whether there should be refinements in
#' the lower and upper tail (see \code{\link{getBreaks}}).
#' @param probs numeric vector with values in \eqn{[0, 1]}; if \code{method} is
#' \code{"multinom"} and \code{breaks} is \code{NULL}, this gives probabilities
#' for quantiles to be used as (positive) break points. If supplied, this is
#' preferred over \code{equidist}.
#' @param gpd logical; if \code{method} is \code{"multinom"}, this indicates
#' whether the upper tail of the personal income should be simulated by random
#' draws from a (truncated) generalized Pareto distribution rather than a
#' uniform distribution.
#' @param threshold a numeric value; if \code{method} is \code{"multinom"},
#' values for categories above \code{threshold} are drawn from a (truncated)
#' generalized Pareto distribution.
#' @param est a character string; if \code{method} is \code{"multinom"}, the
#' estimator to be used to fit the generalized Pareto distribution.
#' @param const numeric; if \code{method} is \code{"twostep"}, this gives a
#' constant to be added before log transformation.
#' @param alpha numeric; if \code{method} is \code{"twostep"}, this gives
#' trimming parameters for the sample data. Trimming is thereby done with
#' respect to the variable specified by \code{additional}. If a numeric vector
#' of length two is supplied, the first element gives the trimming proportion
#' for the lower part and the second element the trimming proportion for the
#' upper part. If a single numeric is supplied, it is used for both. With
#' \code{NULL}, trimming is suppressed.
#' @param residuals logical; if \code{method} is \code{"twostep"}, this
#' indicates whether the random error terms should be obtained by draws from
#' the residuals. If \code{FALSE}, they are drawn from a normal distribution
#' (median and MAD of the residuals are used as parameters).
#' @param components a character vector specifying the income components in
#' \code{dataS} (to be simulated for the population data).
#' @param conditional an optional character vector specifying categorical
#' contitioning variables for resampling of the income components. The
#' fractions occurring in \code{dataS} are then drawn from the respective
#' subsets defined by these variables.
#' @param keep a logical indicating whether variables computed internally in
#' the procedure (such as the original IDs of the corresponding households in
#' the underlying sample, age categories or income categories) should be stored
#' in the resulting population data.
#' @param maxit,MaxNWts control parameters to be passed to
#' \code{\link[nnet]{multinom}} and \code{\link[nnet]{nnet}}. See the help file
#' for \code{\link[nnet]{nnet}}.
#' @param tol if \code{method} is \code{"twostep"}, a small positive numeric
#' value or \code{NULL} (see \code{\link{simContinuous}}).
#' @param seed optional; an integer value to be used as the seed of the random
#' number generator, or an integer vector containing the state of the random
#' number generator to be restored.
#' @return An object of class \code{\linkS4class{simPopObj}} containing the
#' simulated EU-SILC population data as well as the underlying sample.
#' @note This is a wrapper calling \code{\link{simStructure}},
#' \code{\link{simCategorical}}, \code{\link{simContinuous}} and
#' \code{\link{simComponents}}.
#' @author Andreas Alfons and Stefan Kraft and Bernhard Meindl
#' @seealso \code{\link{simStructure}}, \code{\link{simCategorical}},
#' \code{\link{simContinuous}}, \code{\link{simComponents}}
#' @keywords datagen
#' @export
#' @examples
#' 
#' data(eusilcS) # load sample data
#' 
#' \dontrun{
#' ## long computation time
#' # multinomial model with random draws
#' eusilcM <- simEUSILC(eusilcS, upper = 200000, equidist = FALSE)
#' summary(eusilcM)
#' 
#' # two-step regression
#' eusilcT <- simEUSILC(eusilcS, method = "twostep")
#' summary(eusilcT)
#' }
#' 
simEUSILC <- function(dataS, hid = "db030", wh = "db090",
  wp = "rb050", hsize = NULL, strata = "db040",
  pid = NULL, age = "age", gender = "rb090",
  categorizeAge = TRUE, breaksAge = NULL,
  categorical = c("pl030", "pb220a"),
  income = "netIncome", method = c("multinom", "twostep"),
  breaks = NULL, lower = NULL, upper = NULL,
  equidist = TRUE, probs = NULL, gpd = TRUE,
  threshold = NULL, est = "moments", const = NULL,
  alpha = 0.01, residuals = TRUE,
  components = c("py010n", "py050n", "py090n", "py100n", "py110n", "py120n", "py130n", "py140n"),
  conditional = c(getCatName(income), "pl030"),
  keep = TRUE, maxit = 500, MaxNWts = 1500,
  tol = .Machine$double.eps^0.5, seed) {

  ## initializations
  if ( !missing(seed) ) {
    set.seed(seed)
  }
  if ( !is.character(age) || length(age) != 1 ) {
    stop("'age' must be a character string specifying exactly one column of 'dataS'.\n")
  }
  if ( !is.character(gender) || length(gender) != 1 ) {
    stop("'gender' must be a character string specifying exactly one column of 'dataS'.\n")
  }
  if(!is.character(income) || length(income) != 1) {
    stop("'income' must be a character string specifying exactly one column of 'dataS'\n.")
  }

  ## simulate household structure
  structure <- c(age, gender)
  inp <- specifyInput(data=dataS, hhid=hid, weight=wh, hhsize=hsize, strata=strata, pid=pid)
  simPop <- simStructure(dataS=inp, basicHHvars=structure)

  dataS <- simPop@sample@data
  dataP <- simPop@pop@data

  # construct age categories (if requested)
  categorizeAge <- isTRUE(categorizeAge)
  if ( categorizeAge ) {
    ageCat <- getCatName(age)
    # check breakpoints
    if ( is.null(breaksAge) ) {
      r <- c(range(dataS[[age]], na.rm=TRUE))
      s <- seq(15, 80, 5)
      breaksAge <- c(r[1], s[s > r[1] & s < r[2]], r[2])
    } else if ( !is.numeric(breaksAge) || length(breaksAge) < 2 ) {
      stop("'breaksAge' must be a numeric vector with length >= 2")
    }
    # categorize
    dataS[[ageCat]] <- as.character(cut(dataS[[age]], breaks=breaksAge, include.lowest=TRUE))
    dataP[[ageCat]] <- as.character(cut(dataP[[age]], breaks=breaksAge, include.lowest=TRUE))
    # use age class as predictor instead of age
    structure <- c(ageCat, gender)
  } else ageCat <- NULL
  simPop@pop@data <- dataP
  simPop@sample@data <- dataS

  ## simulate additional categorical variables
  basic <- c(structure, if(is.null(hsize)) "hsize" else hsize)
  simPop <- simCategorical(simPop, additional=categorical, maxit=maxit, MaxNWts=MaxNWts)

  dataS <- simPop@sample@data
  dataP <- simPop@pop@data

  ## simulate income
  basic <- union(basic, categorical)
  method <- match.arg(method)
  useMultinom <- method == "multinom"
  zeros <- TRUE
  # compute income categories
  exclude <- getExclude(dataS[, c(wp, strata, basic, income), with=F])
  if ( length(exclude) > 0 ) {
    incomeS <- dataS[[income]][-exclude]
    wpS <- dataS[[wp]][-exclude]
  } else {
    incomeS <- dataS[[income]]
    wpS <- dataS[[wp]]
  }
  missingBreaks <- is.null(breaks)
  if ( missingBreaks ) {
    if ( is.null(upper) && gpd ) {
      upper <- Inf
    }
    breaks <- getBreaks(incomeS, wpS, zeros, lower, upper, equidist, probs)
  }
  incomeCat <- getCatName(income)
  dataS[[incomeCat]] <- getCat(dataS[[income]], breaks, zeros)
  simPop@sample@data <- dataS
  simPop@pop@data <- dataP

  if ( useMultinom ) {
    # multinomial model with random draws from resulting categories
    if ( is.null(threshold) && missingBreaks && (!isTRUE(equidist) || !is.null(probs)) ) {
      threshold <- breaks[length(breaks)-2]
    }
    simPop <- simContinuous(simPop, additional=income, zeros=zeros, breaks=breaks,
      gpd=gpd, threshold=threshold, est=est, keep=TRUE, maxit=maxit, MaxNWts=MaxNWts)
  } else {
    # two-step model
    simPop <- simContinuous(simPop, additional=income, method="lm", zeros=zeros,
      breaks=breaks, log=TRUE, const=const, alpha=alpha,
      residuals=residuals, maxit=maxit, MaxNWts=MaxNWts, tol=tol)
    dataP <- simPop@pop@data
    dataP[[incomeCat]] <- getCat(dataP[[income]], breaks, zeros)
    simPop@pop@data <- dataP
  }

  ## simulate income components
  simPop <- simComponents(simPop, total=income, components=components, conditional=conditional)

  # round income components and adjust income
  dataP <- simPop@pop@data
  dataP[, (components) := lapply(.SD, round, digits=2), .SDcols = components]
  dataP[[income]] <- rowSums(dataP[, components, with=F])

  # return population data
  if ( !keep ) {
    dataP <- dataP[, setdiff(names(dataP), c(ageCat, incomeCat)), with=F]
  }
  simPop@pop@data <- dataP
  return(invisible(simPop))
}

