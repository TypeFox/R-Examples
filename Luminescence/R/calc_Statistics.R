#' Function to calculate statistic measures
#'
#' This function calculates a number of descriptive statistics for De-data,
#' most fundamentally using error-weighted approaches.
#'
#'
#' @param data \code{\link{data.frame}} or \code{\linkS4class{RLum.Results}}
#' object (required): for \code{data.frame} two columns: De (\code{data[,1]})
#' and De error (\code{data[,2]}). To plot several data sets in one plot the
#' data sets must be provided as \code{list}, e.g. \code{list(data.1, data.2)}.
#'
#' @param weight.calc \code{\link{character}}: type of weight calculation. One
#' out of \code{"reciprocal"} (weight is 1/error), \code{"square"} (weight is
#' 1/error^2). Default is \code{"square"}.
#'
#' @param digits \code{\link{integer}} (with default): round numbers to the specified digits. If
#' digits is set to \code{NULL} nothing is rounded.
#'
#' @param na.rm \code{\link{logical}} (with default): indicating whether NA
#' values should be stripped before the computation proceeds.
#'
#' @return Returns a list with weighted and unweighted statistic measures.
#'
#' @section Function version: 0.1.4
#'
#' @author Michael Dietze, GFZ Potsdam (Germany)
#'
#' @examples
#'
#' ## load example data
#' data(ExampleData.DeValues, envir = environment())
#'
#' ## show a rough plot of the data to illustrate the non-normal distribution
#' plot_KDE(ExampleData.DeValues$BT998)
#'
#' ## calculate statistics and show output
#' str(calc_Statistics(ExampleData.DeValues$BT998))
#'
#' ## now the same for 10000 normal distributed random numbers with equal errors
#' x <- as.data.frame(cbind(rnorm(n = 10^5, mean = 0, sd = 1),
#'                          rep(0.001, 10^5)))
#'
#' ## note the congruent results for weighted and unweighted measures
#' str(calc_Statistics(x))
#'
#' @export
calc_Statistics <- function(
  data,
  weight.calc = "square",
  digits = NULL,
  na.rm = TRUE
) {
  ## Check input data
  if(is(data, "RLum.Results") == FALSE &
       is(data, "data.frame") == FALSE) {
    stop(paste("[calc_Statistics()] Input data format is neither",
               "'data.frame' nor 'RLum.Results'"))
  } else {
    if(is(data, "RLum.Results") == TRUE) {
      data <- get_RLum(data)[,1:2]
    }
  }

  ##strip na values
  if(na.rm){
    data <- na.exclude(data)
  }

  ## handle error-free data sets
  if(ncol(data) == 1) {
    data <- cbind(data, rep(NA, length(data)))
  }

  ## replace Na values in error by 0
  data[is.na(data[,2]),2] <- 0

  if(sum(data[,2]) == 0) {
    warning("All errors are NA or zero! Automatically set to 10^-9!")
    data[,2] <- rep(x = 10^-9, length(data[,2]))
  }

  if(weight.calc == "reciprocal") {
    S.weights <- 1 / data[,2]
  } else if(weight.calc == "square") {
    S.weights <- 1 / data[,2]^2
  } else {
    stop ("[calc_Statistics()] Weight calculation type not supported!")
  }

  S.weights <- S.weights / sum(S.weights)

  ## calculate n
  S.n <- nrow(data)

  ## calculate mean
  S.mean <- mean(x = data[,1],
                 na.rm = na.rm)

  S.wg.mean <- weighted.mean(x = data[,1],
                             w = S.weights,
                             n.rm = na.rm)

  ## calculate median
  S.median <- median(x = data[,1],
                     na.rm = na.rm)

  S.wg.median <- S.median

  ## calculate absolute standard deviation
  S.sd.abs <- sd(x = data[,1],
                 na.rm = na.rm)

  S.wg.sd.abs <- sqrt(sum(S.weights * (data[,1] - S.wg.mean)^2) /
                        (((S.n - 1) * sum(S.weights)) / S.n))

  ## calculate relative standard deviation
  S.sd.rel <- S.sd.abs / S.mean * 100

  S.wg.sd.rel <- S.wg.sd.abs / S.wg.mean * 100

  ## calculate absolute standard error of the mean
  S.se.abs <- S.sd.abs / sqrt(S.n)

  S.wg.se.abs <- S.wg.sd.abs / sqrt(S.n)

  ## calculate relative standard error of the mean
  S.se.rel <- S.se.abs / S.mean * 100

  S.wg.se.rel <- S.wg.se.abs / S.wg.mean * 100

  ## calculate skewness
  S.skewness <- 1 / S.n * sum(((data[,1] - S.mean) / S.sd.abs)^3)

  ## calculate kurtosis
  S.kurtosis <- 1 / S.n * sum(((data[,1] - S.mean) / S.sd.abs)^4)

  ## create list objects of calculation output
  S.weighted <- list(n = S.n,
                     mean = S.wg.mean,
                     median = S.wg.median,
                     sd.abs = S.wg.sd.abs,
                     sd.rel = S.wg.sd.rel,
                     se.abs = S.wg.se.abs,
                     se.rel = S.wg.se.rel)


  if(!is.null(digits)){
     S.weighted <- sapply(names(S.weighted), simplify = FALSE, USE.NAMES = TRUE, function(x){
      round(S.weighted[[x]], digits = digits)

    })

  }

  S.unweighted <- list(n = S.n,
                       mean = S.mean,
                       median = S.median,
                       sd.abs = S.sd.abs,
                       sd.rel = S.sd.rel,
                       se.abs = S.se.abs,
                       se.rel = S.se.rel,
                       skewness = S.skewness,
                       kurtosis = S.kurtosis)

  if(!is.null(digits)){
    S.unweighted  <- sapply(names(S.unweighted), simplify = FALSE, USE.NAMES = TRUE, function(x){
      round(S.unweighted [[x]], digits = digits)

    })

  }

  list(weighted = S.weighted,
       unweighted = S.unweighted)

}

