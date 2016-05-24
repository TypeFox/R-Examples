#' Checks if an object can be used to perform reproduction toxicity data analysis
#'
#' The \code{reproDataCheck} function can be used to check if an object
#' containing data from a reproduction toxicity assay meets the expectations
#' of the function \code{\link{reproData}}.
#'
#' Since in MORSE reproduction datasets are a special case of survival datasets,
#' \code{reproDataCheck} performs the same verifications than
#' \code{\link{survDataCheck}} plus additional ones that are specific to
#' reproduction data.
#'
#' @aliases reproDataCheck
#'
#' @param data any object
#' @param diagnosis.plot if \code{TRUE}, may produce a diagnosis plot
#'
#' @return The function returns a \code{data.frame} similar to the one returned
#' by \code{\link{survDataCheck}}, except that it may contain the following
#' additional error \code{id}s:
#' \itemize{
#' \item \code{NreproInteger}: column \code{Nrepro} contains values of class other than \code{integer}
#' \item \code{Nrepro0T0}: \code{Nrepro} is not 0 at time 0 for each concentration and each replicate
#' \item \code{Nsurvt0Nreprotp1P}: at a giving time \eqn{T}, the number of
#' alive individuals is null and the number of collected offspring is not null
#' for the same replicate and the same concentration at time \eqn{T+1}
#' }
#'
#' @note If an error of type \code{dataframeExpected} or \code{missingColumn}
#' is detected, the function
#' \code{reproDataCheck} is stopped. When no error is detected the
#' \code{reproDataCheck}
#' function returns an empty dataframe.
#'
#' @seealso \code{\link{reproData}}
#'
#' @keywords check
#'
#' @examples
#'
#' # Run the check data function
#' data(copper)
#' reproDataCheck(copper)
#'
#' # Now we insert an error in the dataset, by setting a non-zero number of
#' # offspring at some time, although there is no surviving individual in the
#' # replicate from the previous time point.
#' copper$Nrepro[148] <- as.integer(1)
#' reproDataCheck(copper)
#'
#' @export
reproDataCheck <- function(data, diagnosis.plot = TRUE) {

  ##
  ## 1 run the tests of survDataCheck
  ##
  errors <- survDataCheck(data, diagnosis.plot = FALSE)
  if ("dataframeExpected" %in% errors$id || "missingColumn" %in% errors$id)
    return(errors)

  ##
  ## 1' test if the column names "Nrepro" exists
  ##
  if (! "Nrepro" %in% colnames(data)) {
    msg <- "The column Nrepro is missing."
    errors <- errorTableAdd(errors, "missingColumn", msg)
    return(errors)
  }

  ##
  ## 2' test if Nrepro is integer
  ##
  if (!is.integer(data$Nrepro)) {
    msg <- "Column 'Nrepro' must contain only integer values."
    errors <- errorTableAdd(errors, "NreproInteger", msg)
  }

  ##
  ## 3'test Nrepro = 0 at time 0
  ##
  datatime0 <- data[data$time == 0, ] # select data for initial time points
  if (any(datatime0$Nrepro > 0)) { # test if Nrepro > 0 at time 0
    msg <- "Nrepro should be 0 at time 0 for each concentration and each replicate."
    errors <- errorTableAdd(errors, "Nrepro0T0", msg)
  }

  subdata <- split(data, list(data$replicate, data$conc), drop = TRUE)
  consistency <- function(subdata) {
    # Function to be used on a subdataset corresponding to one replicate at one
    # concentration.
    # This function checks:
    #   - if at each time T for which Nsurv = 0, Nrepro = 0 at time T+1

    # errors consitency dataframe
    errors <- errorTableCreate()

    ##
    ## 4' test Nsurv = 0 at time t => Nrepro > 0 at time t-1
    ##
    NsurvT <- subdata$Nsurv[-length(subdata$Nsurv)]
    NreproTplus1 <- subdata$Nrepro[-1]
    if (any(NreproTplus1[NsurvT == 0] > 0)) {
      msg <- paste("For replicate ",
                   unique(subdata$replicate),
                   " and concentration ", unique(subdata$conc),
                   ", there are some Nsurv = 0 followed by Nrepro > 0 at the next time point.",
                   sep = "")
      errors <- errorTableAdd(errors, "Nsurvt0Nreprotp1P", msg)
    }
    return(errors)
  }
  res <- by(data, list(data$replicate, data$conc), consistency)
  consistency.errors <- do.call("errorTableAppend", res)
  errors <- errorTableAppend(errors, consistency.errors)

  if (diagnosis.plot && "NsurvIncrease" %in% errors$id) {
    survDataPlotFull(data, ylab = "Number of surviving individuals")
      }

  return(errors)
}
