##' Response and correlation function analysis
##' 
##' This function calculates (potentially moving) response and
##' correlation functions from tree-ring chronologies and monthly
##' climatic data. For the moving case, the calculation is performed
##' repeatedly for consecutive time windows. Function parameters may
##' be bootstrapped to calculate their significance and confidence
##' intervals.
##' @details This function builds upon and extents the functionality
##' of programme DENDROCLIM2002 (Biondi and Waikul, 2004), and will
##' calculate bootstrapped (and non-bootstrapped) moving and static
##' response and correlation functions in a similar fashion as
##' described in the above mentioned paper. Important extensions
##' include a very flexible parameter selection model (see below), the
##' possibility to use an unlimited number of climate parameters, and
##' the option to use exact bootstrapping.
##'   
##' Input chronology data can be a \code{data.frame} such as produced
##' by function \code{chron} of package dplR. It has to be a
##' \code{data.frame} with at least one column containing the tree-ring
##' indices, and the corresponding years as \code{rownames}.
##'   
##' For climatic input data, there are three possibilities: Firstly,
##' input climatic data can be a \code{data.frame} or \code{matrix}
##' consisting of at least 3 rows for years, months and at least one
##' climate parameter in the given order. Secondly, input climatic
##' data can be a single \code{data.frame} or \code{matrix} in the
##' style of the original DENDROCLIM2002 input data, i.e. one
##' parameter with 12 months in one row, where the first column
##' represents the year. Or thirdly, input climatic data can be a list
##' of several of the latter described \code{data.frame} or
##' \code{matrices}. As an internal format dispatcher checks the
##' format automatically, it is absolutely necessary that in all three
##' cases, only complete years (months 1-12) are provided. It is not
##' possible to mix different formats in one go.
##'   
##' Parameters can be selected with the 'selection' parameter in two
##' different ways. The default value is -6:9. This is equivalent to
##' the standard settings in DENDROCLIM2002 and bootRes, and selects
##' from all climate variables all months from previous year's June
##' (-6, previous year's months are specified as negative integers) to
##' current years September (9, months of the current year are
##' specified as positive integers) as model parameters.
##'   
##' More complex parameter selections can be obtained by the
##' \emph{modifiers} provided in treeclim: \code{.range},
##' \code{.mean}, and \code{.sum}.  \code{.range} corresponds the
##' example above, where all specified months are used, while
##' \code{.sum} and \code{.mean} will use the sums and means of the
##' specified months. These modifiers also allow to select specific
##' climatic variables, addressed by name. Thus, \code{.mean(4:8,
##' "temp")} will select the mean for climate parameter "temp" for the
##' months April to August. Not only ranges, but also individual
##' vectors can be used for month specification, like e.g.,
##' \code{.range(c(1, 3, 4, 5)}.
##'   
##' The modifiers can be chained together using the '+' symbol, which
##' makes it possible to create arbitrarily complex selections of
##' climate parameters for calibration.  E.g., \code{.mean(2:5,
##' "temp") + .sum(2:5, "prec")} will yield the February-to-May mean
##' for the variable "temp" and the sum of the variable "prec" for the
##' same time. While there is no limitation for number of lists that
##' can be chained together, 'dcc' will not check for meaningful
##' specifications. Testing smart hypotheses is up the researcher.
##'   
##' For the exclusion of months, the convenience function
##' \code{excludefrom()} (or short \code{exfr()}) is provided. E.g.,
##' \code{.range(excludefrom(-6:10, -11:3))} will yield the monthly
##' values of all parameters for the months previous June (-6) to
##' current October (10), but without the months previous November
##' (-11) to current March (3) in between. While it is also possible
##' to supply arbitrary vectors as month specification, and not only
##' ranges as shown in most of the examples here, this way of
##' excluding e.g. the dormant season is far more convenient.
##'   
##' 1000 bootstrap samples are taken from the original distributions
##' of climate and tree-ring data, either using the stationary
##' bootstrap (Politis and Romano 1994, \code{boot = "stationary"}) or
##' classical bootstrap (DENDROCLIM2002-style, \code{boot =
##' "std"}). The stationary bootstrap mimics the stationary properties
##' of the original time series in the resampled time series by
##' resampling within blocks. Within each block, the number of
##' observations is random and has a geometric
##' distribution. Consequently, the choice of the distribution
##' parameter will affect the autocorrelation structure of the
##' resampled time series. Optimal (expected) block length is chosen
##' according to Politis and White (2004). In the case of response
##' function analysis, an eigen decomposition of the standardized
##' predictor matrix is performed. Nonrelevant eigenvectors are
##' removed using the PVP criterion (Guiot, 1990), principal component
##' scores are then calculated from the matrices of reduced
##' eigenvectors and standardized climatic predictors. Response
##' coefficients are found via singular value decomposition, and
##' tested for significance using the 95\% percentile range method
##' (Dixon, 2001). In case of correlation function analysis, the
##' coefficients are Pearson's correlation coefficients. The same
##' method for significance testing is applied.
##'   
##' There is also the option to use exact bootstrapping like
##' implemented in seascorr (Meko et al. 2011, \code{boot =
##' "exact"}). In this case, circulant embedding is used to simulate
##' the tree-ring data 1000 times as time series with the same
##' frequency characteristics like the original time-series (Percival
##' & Constantine, 2006). Empirical non-exceedence probabilities are
##' used to test the coefficients of the response/correlation function
##' with the original data for significance. For the exact
##' bootstrapping case, no confidence intervals for the
##' response/correlation coefficients can be computed.
##' @param chrono \code{data.frame} containing a tree-ring
##' chronologies, e.g. as obtained by \code{chron} of package dplR.
##' @param climate either a \code{data.frame} or \code{matrix} with
##' climatic data in monthly resolution, with year, month and climate
##' parameters in columns (all columns except year and month will be
##' recognized as parameters for response or correlation function), or
##' a single \code{data.frame} or \code{matrix} in 13-column format
##' (see below), or list of several of the latter.
##' @param selection either a numeric vector, a modifier, or a chain
##' of modifiers specifying the parameter selection for the model (see
##' Details).
##' @param method \code{character} string specifying the calculation
##' method.  Possible values are \dQuote{response} and
##' \dQuote{correlation}. Partial strings are ok.
##' @param moving \code{logical}; should the analyis be carried out in
##' moving windows. Defaults to \code{FALSE}.
##' @param win_size integer giving the window size for each
##' recalculation in years.
##' @param win_offset integer giving the number of years between each
##' window start in years.
##' @param start_last \code{logical} flag indicating whether the first
##' window should start at the rear end (youngest part of the series)
##' or not.
##' @param timespan \code{integer} vector of length 2 specifying the
##' time interval (in years) to be considered for analysis. Defaults
##' to the maximum possible interval.
##' @param var_names \code{character} vector with variable
##' names. Defaults to corresponding column names of \code{data.frame}
##' clim.
##' @param ci \code{numerical} value to set the test level for
##' significance test (values 0.01, 0.05 and 0.1 are allowed); the
##' confidence intervals are adapted accordingly.
##' @param boot \code{character} indicating which bootstrap method
##' should be used, one of \code{c("stationary", "std", "exact")}
##' @param sb \code{logical} flag indicating whether textual status
##' bar for moving case should be suppressed. Suppression is
##' recommended for e.g.  Sweave files.
##' @return 'dcc' returns an 'object' of class '"tc_dcc"'.
##'
##' The functions 'summary' and 'plot' are used to obtain and print a
##' summary of the results, and to create a plot. The function 'coef'
##' can be used to extract the coefficients.
##'
##' An object of class '"tc_dcc"' is a list containing at least the
##' following components:
##'
##' \item{call}{the call made to function 'dcc'}
##' 
##' \item{coef}{the coefficients, themselves being an object of class
##' 'tc_coef' for the static case, and of class 'tc_mcoef' for the
##' moving case. Objects of class 'tc_coef' are single data.frames,
##' while objects of class 'tc_mcoef' are lists of seperate
##' data.frames for the coefficients ($coef), upper and lower
##' confidence interval ($ci_upper and $ci_lower), and significance
##' flags ($significant)}
##' 
##' \item{design}{the design matrix on which this call to 'dcc'
##' operates}
##' 
##' \item{truncated}{the input data truncated to the common timespan
##' or the specified timespan}
##' 
##' \item{original}{the original input data, with the climate data
##' being recast into a single data.frame}
##' @references Biondi, F & Waikul, K (2004) DENDROCLIM2002: A C++
##' program for statistical calibration of climate signals in
##' tree-ring chronologies.  \emph{Computers & Geosciences} 30:303-311
##'   
##' Dixon, PM (2001) Bootstrap resampling. In: El-Shaarawi, AH,
##' Piegorsch, WW (Eds.), \emph{The Encyclopedia of
##' Environmetrics}. Wiley, New York.
##'   
##' Guiot, J (1991) The boostrapped response function. \emph{Tree-Ring
##' Bulletin} 51:39-41
##'   
##' Meko DM, Touchan R, Anchukaitis KJ (2011) Seascorr: A MATLAB
##' program for identifying the seasonal climate signal in an annual
##' tree-ring time series.  \emph{Computers \& Geosciences}
##' 37:1234-241
##'   
##' Percival DB, Constantine WLB (2006) Exact simulation of Gaussian
##' Time Series from Nonparametric Spectral Estimates with Application
##' to Bootstrapping. \emph{Statistics and Computing} 16:25-35
##'
##' Patton, A. and D.N. Politis and H. White (2009), "CORRECTION TO
##' Automatic block-length selection for the dependent bootstrap" by
##' D. Politis and H. White", Econometric Reviews 28(4), 372-375.
##' 
##' Politis, D.N. and H. White (2004), Automatic block-length
##' selection for the dependent bootstrap, Econometric Reviews 23(1),
##' 53-70.
##' @examples
##' \dontrun{
##' dc_resp <- dcc(muc_spruce, muc_clim)
##' }
##' @author Christian Zang; the original MATLAB code for exact
##' bootstrapping was written by Dave Meko
##' @export
dcc <- function(chrono,
               climate,
               selection = -6:9,
               method = "response",
               moving = FALSE,
               win_size = 25,
               win_offset = 1,
               start_last = TRUE,
               timespan = NULL,
               var_names = NULL,
               ci = 0.05,
               boot = "stationary",
               sb = TRUE
               )
{

  ## climate data are correctly formatted
  climate <- as_tcclimate(climate)
  ## when var_names are supplied, apply appropriately
  if (!is.null(var_names)) {
    varno <- dim(climate)[2] - 2
    if (length(var_names) != varno) {
      stop("Count of supplied variable names does not match count of variables in climate data.")
    } else {
      names(climate)[3:(dim(climate)[2])] <- var_names
    }
  }

  ## check the ci
  if (!any(ci == c(0.1, 0.05, 0.01)))
    stop("'ci' has to be one of [0.1, 0.05, 0.01].")


  ## check the month selection
  ## when numeric vector = "classic" (bR 1.X) parameter selection
  if (is.numeric(selection)) {
    selection <- eval(
      substitute(
        .range(.months = sel, .variables = NULL),
        list(sel = selection)
        )
      )
  }
  
  ## check, if we have months in the correct numerical representation, and
  ## preevaluate month specification for correct data truncation
  monthcheck <- check_months(selection)
  if (monthcheck$check == FALSE) {
    stop("Please specify months with numbers from -12 (previous december) to 12 (current december).")
  }
  minmonth <- monthcheck$minmonth

  ## truncate climate and tree-ring data to common or specified
  ## time span
  truncated_input <- truncate_input(chrono, climate,
                                   timespan = timespan, minmonth,
                                   moving)

  ## check if the timespan matches with win_size
  if (moving) {
    time_length <- dim(truncated_input$climate)[1]
    if (time_length <= win_size) {
      stop(paste("timespan is shorter than win_size. Consider adapting timespan to at least ",
                 win_size, ", or win_size to less than ", time_length,
                 ".", sep = ""))
    }
  }

  ## generate matrix of climate data for further processing
  climate_pmat <- make_pmat(truncated_input$climate, truncated_input$pad)

  ## generate design matrix for calibration by evaluating the selections
  design <- tc_design(selection, climate_pmat)

  ## check if number of parameters is smaller than number of observations
  n_params <- dim(design$aggregate)[2]
  if (moving) {
    n_obs <- win_size
  } else {
    n_obs <- dim(design$aggregate)[1]  
  }
  if (n_params > n_obs) {
    stop(
      paste("Overlapping time span of chrono and climate records is smaller than number of parameters! Consider adapting the number of parameters to a maximum of ",
            n_obs, ".", sep = ""))
  }

  ## pass chrono and design matrix to the respective analysis functions

  .method <- match.arg(method, c("response", "correlation"))
  .boot <- match.arg(boot, c("stationary", "std", "exact"))
  
  ## static functions

  if (!moving) {
    if (.method == "response") {
      dc <- tc_response(truncated_input$chrono, design,
                       ci = ci,
                       boot = .boot)
    }
    
    if (.method == "correlation") {
      dc <- tc_correlation(truncated_input$chrono, design,
                          ci = ci,
                          boot = .boot)
    }
  }
  
  ## moving functions

  if (moving) {
    dc <- tc_mfunc(truncated_input$chrono, design,
                  ci = ci,
                  sb = sb,
                  method = .method,
                  start_last = start_last,
                  win_size = win_size,
                  win_offset = win_offset,
                  boot = .boot)
  }

  ## return everything in a comprehensible manner
  
  dcc_out <- list()
  dcc_out$call <- match.call()
  dcc_out$coef <- dc$result
  dcc_out$ac <- dc$ac
  dcc_out$design <- design
  dcc_out$truncated <- list(tree = truncated_input$chrono,
                           climate = truncated_input$climate)
  dcc_out$original <- list(tree = chrono,
                          climate = climate)

  class(dcc_out) <- c("tc_dcc", "list")

  dcc_out
}
