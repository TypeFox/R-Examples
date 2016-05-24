##' Seasonal (partial) correlation analysis
##' 
##' Calculate seasonal correlation with primary and secondary climate variables
##' and tree-ring data, similar to the seascorr function for MATLAB.
##' @details This function mimicks the behaviour of the MATLAB function seascorr
##'   (Meko et al. 2011), which calculates partial correlations of tree-ring
##'   data with a primary and a secondary climatic variable for seasons of
##'   different lengths.
##'   
##'   Input chronology data can be a \code{data.frame} such as produced by
##'   function \code{chron} of package dplR. It has to be a \code{data.frame}
##'   with at least one column containing the tree-ring indices, and the
##'   corresponding years as \code{rownames}.
##'   
##'   For climatic input data, there are three possibilities: Firstly, input
##'   climatic data can be a \code{data.frame} or \code{matrix} consisting of at
##'   least 3 rows for years, months and at least one climate parameter in the
##'   given order. Secondly, input climatic data can be a single
##'   \code{data.frame} or \code{matrix} in the style of the original
##'   DENDROCLIM2002 input data, i.e. one parameter with 12 months in one row,
##'   where the first column represents the year. Or thirdly, input climatic
##'   data can be a list of several of the latter described \code{data.frame} or
##'   \code{matrices}. As an internal format dispatcher checks the format
##'   automatically, it is absolutely necessary that in all three cases, only
##'   complete years (months 1-12) are provided. It is not possible to mix
##'   different formats in one go.
##'   
##'   The `complete` parameter specifies the months of the current year in which
##'   tree-growth is assumed to finish. This month marks the last month of the
##'   first season, and starting from here, 14 different seasons are computed
##'   for each specified season length in one-month steps. E.g., for a starting
##'   value of 9 (current September) and season length of 3 months, the first
##'   season comprises current July to current September, the second season 
##'   comprises current June to current August, and the last season comprises
##'   previous June to previous August. This results in 14 seasons for a given
##'   season length. An arbitrary number of season lengths can be specified.
##'   
##'   The choice for primary vs. secondary variable can be made either via
##'   numeric selection (the integer value 1 stands for the first variable in
##'   the supplied climate data set), or by name ("temp", when one of the
##'   variables is named "temp"). The correlation of the primary variable with
##'   tree-growth is computed as the simple (Pearson) correlation coefficient,
##'   while the influence of the secondary variable on tree-growth is computed
##'   with the influence of the primary variable on tree-growth removed.
##'   
##'   Like in the original seascorr program, the significance of each (partial)
##'   correlation is evaluated using exact bootstrapping by circulant embedding
##'   of the tree-ring data (Percival \& Constantine, 2006).
##' 
##' @param chrono \code{data.frame} containing a tree-ring
##'   chronologies, e.g. as obtained by \code{chron} of package dplR.
##' 
##' @param climate either a \code{data.frame} or \code{matrix} with
##'   climatic data in monthly resolution, with year, month and
##'   climate parameters in columns (all columns except year and month
##'   will be recognized as parameters for response or correlation
##'   function), or a single \code{data.frame} or \code{matrix} in
##'   13-column format (see below), or a list of several of the
##'   latter.
##' 
##' @param var_names \code{character} vector with variable
##'   names. Defaults to corresponding column names of
##'   \code{data.frame} clim.
##' 
##' @param timespan \code{integer} vector of length 2 specifying the
##'   time interval (in years) to be considered for analysis. Defaults
##'   to the maximum possible interval.
##'
##' @param complete \code{integer} scalar, month when tree-ring growth
##'   is expected to have finished.
##'
##' @param season_lengths \code{numeric} vector giving the lengths of
##'   the seasons for variable grouping
##'
##' @param primary position \code{numeric} or name \code{character} of
##'   primary climate variable
##' 
##' @param secondary position \code{numeric} or name \code{character}
##'   of secondary climate variable
##' 
##' @param ci \code{numeric} scalar to set the test level for
##'   significance test (values 0.01, 0.05 and 0.1 are allowed); the
##'   confidence intervals are adapted accordingly.
##' 
##' @return 'seascorr' returns an 'object' of class '"tc_seascorr"'.
##'
##' The 'plot' function is used to obtain a plot of the results.
##'
##' An object of class '"tc_seascorr"' is a list containing at least
##'   the following components:
##'
##' \item{call}{the call made to 'seascorr'}
##' 
##' \item{seasons}{a list of length n, where n is the number of season
##'   lengths provided; each list element consists of a data.frame
##'   with end month, correlation coefficient and significance flag}
##' 
##' \item{truncated}{the input data truncated to the common timespan
##'   or the specified timespan}
##' 
##' \item{original}{the original input data, with the climate data
##'   being recast into a single data.frame}
##'
##' @references
##'
##' Meko DM, Touchan R, Anchukaitis KJ (2011) Seascorr: A MATLAB
##'   program for identifying the seasonal climate signal in an annual
##'   tree-ring time series. Computers & Geosciences, 37, 1234-1241.
##'
##' Percival DB, Constantine WLB (2006) Exact simulation of Gaussian
##'   Time Series from Nonparametric Spectral Estimates with
##'   Application to Bootstrapping. Statistics and Computing 16:25-35
##'
##' @examples
##' sc <- seascorr(muc_fake, muc_clim)
##' sc
##' plot(sc)
##' @author Christian Zang; the procedure incl. exact bootstrapping
##' was was implemented first by Dave Meko in MA
##' @export
seascorr <- function(chrono, climate, var_names = NULL, timespan =
                     NULL, complete = 9, season_lengths = c(1, 3, 6),
                     primary = 1, secondary = 2, ci = 0.05) {

  ## check input
  if (!any(1:12 == complete))
    stop("`complete` must be an integer value between 1 and 12.")

  if (!all(sapply(season_lengths, function(x) any(1:12 == x))))
    stop("`season_lengths` must be a vector of integers between 1 and 12.")

  if (!any(c(0.01, 0.05, 0.1) == ci))
    stop("`ci` must be any of 0.01, 0.05, or 0.1.")
  
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

  ## Makes no sense for less than 2 climate variables
  if (dim(climate)[2] < 4) {
    stop("Two climate variables needed for season correlations.")
  }

  ## get names of primary and secondary variable
  if (any(1:100 == primary)) {
    ## specified by position
    pos <- primary + 2
    if (pos > length(names(climate))) {
      stop("Position for primary variable does not exist.")
    }
    primary_name <- names(climate)[pos]
  } else {
    if (is.character(primary)) {
      ## specified by name
      if (any(names(climate) == primary)) {
        primary_name <- primary
      } else {
        stop("Name for primary variable does not exist.")
      }
    } else {
      stop("Please specify either position or name for primary variable.")
    }
  }

  if (any(1:100 == secondary)) {
    ## specified by position
    pos <- secondary + 2
    if (pos > length(names(climate))) {
      stop("Position for secondary variable does not exist.")
    }
    secondary_name <- names(climate)[pos]
  } else {
    if (is.character(secondary)) {
      ## specified by name
      if (any(names(climate) == secondary)) {
        secondary_name <- secondary
      } else {
        stop("Name for secondary variable does not exist.")
      }
    } else {
      stop("Please specify either position or name for secondary variable.")
    }
  }

  ## check if identical
  if (primary_name == secondary_name) {
    stop("Primary and secondary variable are identical.")
  }

  ## truncate climate and tree-ring data to common or specified
  ## time span
  truncated_input <- truncate_input(chrono, climate,
                                    timespan = timespan, 1,
                                    moving = FALSE)

  ## create raw parameter matrix
  if (truncated_input$pad) {
    # we have no climate data for previous year, but use maximum overlap -> cut
    # tree data accordingly
    truncated_input$chrono <- truncated_input$chrono[-1]
  }
  
  m <- length(truncated_input$chrono)
  
  ## stop if less than 31 years (needed for exact bootstrapping to
  ## work correctly)
  if (m < 31)
    stop("Seasonal correlation analysis needs at least 31 years of data overlap for proxy and climate data.")
  
  pmat <- make_pmat(truncated_input$climate, pad = FALSE)

  ## create seasons, a list entry for each season_length
  seasons1 <- seasons2 <- list()
  n <- length(season_lengths)
  first_month <- -1 * complete + 1
  last_month <- complete

  ## check if largest season spec is feasible
  if (first_month + max(season_lengths) > 0) {
    maxval <- abs(first_month)
    stop(paste("Largest season length is too long. Maximum value is ",
               maxval, ".", sep = ""))
  }

  lmonths <- c(-1:-12, 1:12)
  last_month_index <- which(lmonths == last_month)

  for (i in 1:n) {
    .season_length <- season_lengths[i]
    seasons1[[i]] <- matrix(NA, ncol = 14, nrow = m)
    seasons2[[i]] <- matrix(NA, ncol = 14, nrow = m)

    for (j in 1:14) {
      end <- last_month_index + 1 - j
      start <- end - .season_length + 1
      end_m <- lmonths[end]
      start_m <- lmonths[start]

      ## create parameter selection list for use with eval_selection
      if (start_m == end_m) {
        method <- "full"
      } else {
        method <- "mean"
      }
      primary_selection <- list(method, start_m:end_m, primary_name)
      secondary_selection <- list(method, start_m:end_m, secondary_name)

      seasons1[[i]][,j] <-
        eval_selection(pmat, primary_selection)$aggregate[,1]
      seasons2[[i]][,j] <-
        eval_selection(pmat, secondary_selection)$aggregate[,1]
    }
  }

  chrono_boot <- init_boot_data(pmat, truncated_input$chrono, 1000,
                                "exact")$chrono

  results <- list()
  results$coef <- list()

  for (i in 1:n) {

    results$coef[[i]] <- list()
    
    params <- .Call("treeclim_pcor", PACKAGE = 'treeclim',
                    seasons1[[i]], seasons2[[i]],
                    chrono_boot, truncated_input$chrono)

    results$coef[[i]]$primary <- ptest(params$primary[,2:1001], ci,
                                       params$primary[,1],
                                       "weibull")[,1:2]
    results$coef[[i]]$secondary <- ptest(params$exact[,2:1001], ci,
                                         params$secondary[,1],
                                         "weibull")[,1:2]
  }

  results$call <- match.call()
  results$seasons <- list(
    primary = seasons1,
    secondary = seasons2
    )
  results$truncated <- list(tree = truncated_input$chrono,
                            climate = truncated_input$climate)
  results$original <- list(tree = chrono,
                           climate = climate)
  class(results) <- c("tc_seascorr", "list")
  results
}
