#' Checks if an object can be used to perform survival analysis
#'
#' The \code{survDataCheck} function can be used to check if an object
#' containing survival data is formatted according to the expectations of the
#' \code{survData} function.
#'
#'
#' @aliases survDataCheck
#'
#' @param data any object
#' @param diagnosis.plot if \code{TRUE}, the function may produce diagnosis plots
#'
#' @return The function returns a dataframe with two columns \code{id} and \code{msg} of
#' character strings. When no error is detected this dataframe is empty.
#' Here is the list of possible error \code{id}s and their signification:
#' \tabular{rl}{
#' \code{dataframeExpected} \tab an object of class \code{data.frame} is expected \cr
#' \code{missingColumn} \tab at least one expected column heading is missing \cr
#' \code{firstTime0} \tab the first time point for some (concentration, replicate) is not 0 \cr
#' \code{concNumeric} \tab column \code{conc} contains a value of class other than \code{numeric} \cr
#' \code{timeNumeric} \tab column \code{time} contains a value of class other than \code{numeric} \cr
#' \code{NsurvInteger} \tab column \code{Nsurv} contains a value of class other than \code{integer} \cr
#' \code{tablePositive} \tab some data are negative \cr
#' \code{Nsurv0T0} \tab \code{Nsurv} is 0 at time 0 for some (concentration, replicate) \cr
#' \code{duplicateID} \tab there are two identical (\code{replicate}, \code{conc}, \code{time}) triplets \cr
#' \code{missingReplicate} \tab a replicate is missing at some (time point, concentration) \cr
#' \code{NsurvIncrease} \tab \code{Nsurv} increases at some time point of some (concentration, replicate) \cr
#' \code{ReplicateLabel} \tab replicate labels differ between two time points at some concentration \cr
#' }
#'
#' @note If an error of type \code{dataframeExpected} or \code{missingColumn} is
#' detected, the function \code{survDataCheck} is stopped before looking for
#' other errors.
#'
#' @seealso \code{\link{survData}}
#'
#' @examples
#' # Run the check data function
#' data(zinc)
#' survDataCheck(zinc)
#'
#' # Now we insert an error in the dataset, by artificially increasing the
#' # number of survivors at some time point, in such a way that the number
#' # of indivuals increases in some replicate
#' zinc[25, "Nsurv"] <- as.integer(20)
#' survDataCheck(zinc, diagnosis.plot = TRUE)
#'
#' @importFrom stringr str_c
#' 
#' @export
survDataCheck <- function(data, diagnosis.plot = TRUE) {

  ##
  ## 0. check we have a data.frame
  ##
  if (class(data) != "data.frame") {
    return(errorTableSingleton("dataframeExpected",
                                "A dataframe is expected."))
  }

  ##
  ## 1. assert column names are correct
  ##
  ref.names <- c("replicate","conc","time","Nsurv")
  missing.names <- ref.names[which(is.na(match(ref.names, names(data))))]
  if (length(missing.names) != 0) {
    msg <- paste("The column ", missing.names,
                 " is missing.", sep = "")
    return(errorTableSingleton("missingColumn",msg))
  }

  # Next errors do not prevent from checking others
  errors <- errorTableCreate()

  ##
  ## 2. assert the first time point is zero for each (replicate, concentration)
  ##
  subdata <- split(data, list(data$replicate, data$conc), drop = TRUE)
  if (any(unlist(lapply(subdata, function(x) x$time[1] != 0)))) {
    msg <- "Data are required at time 0 for each concentration and each replicate."
    errors <- errorTableAdd(errors, "firstTime0", msg)
  }

  ##
  ## 3. assert concentrations are numeric
  ##
  if (!is.double(data$conc) && !is.integer(data$conc)) {
    msg <- "Column 'conc' must contain only numerical values."
    errors <- errorTableAdd(errors, "concNumeric", msg)
  }

  ##
  ## 4. assert time is numeric
  ##
  if (!is.numeric(data$time)) {
    msg <- "Column 'time' must contain only numerical values."
    errors <- errorTableAdd(errors, "timeNumeric", msg)
  }

  ##
  ## 5. assert Nsurv contains integer
  ##
  if (!is.integer(data$Nsurv)) {
    msg <- "Column 'Nsurv' must contain only integer values."
    errors <- errorTableAdd(errors, "NsurvInteger", msg)
  }

  ##
  ## 6. assert all data are positive
  ##
  table <- subset(data, select = -c(replicate)) # remove replicate column
  if (any(table < 0.0)) {
    msg <- "Data must contain only positive values."
    errors <- errorTableAdd(errors, "tablePositive", msg)
  }

  ##
  ## 7. assert Nsurv != 0 at time 0
  ##
  datatime0 <- data[data$time == 0, ]  # select data for initial time points
  if (any(datatime0$Nsurv == 0)) { # test if Nsurv != 0 at time 0
    msg <- "Nsurv should be different to 0 at time 0 for each concentration and each replicate."
    errors <- errorTableAdd(errors, "Nsurv0T0", msg)
  }

  ##
  ## 8 assert each (replicate, concentration, time) triplet is unique
  ##
  ID <- idCreate(data) # ID vector
  if (any(duplicated(ID))) {
    msg <- paste("The (replicate, conc, time) triplet ",
                 ID[duplicated(ID)],
                 " is duplicated.", sep = "")
    errors <- errorTableAdd(errors, "duplicatedID", msg)
  }
  consistency <- function(subdata) {
    # Function to be used on a subdataset corresponding to one replicate at one
    # concentration.
    # This function checks:
    #   - if each replicate appears once and only once at each time
    #   - if Nsurv is never increasing with time

    errors <- errorTableCreate()

    ##
    ## 9. assert there is the same number of replicates for each conc and time
    ##
    if (length(subdata$replicate) != length(unique(data$time))) {
      msg <- paste("Replicate ", unique(subdata$replicate),
                   " is missing for at least one time points at concentration ",
                   unique(subdata$conc), ".", sep = "")
      errors <- errorTableAdd(errors, "missingReplicate", msg)
    }

    ##
    ## 10. assert Nsurv never increases with time
    ##
    nsurv.increase <- subdata$Nsurv[-length(subdata$Nsurv)] < subdata$Nsurv[-1]
    if (any(nsurv.increase)) {
      msg <- paste("For replicate ", unique(subdata$replicate),
                   " and concentration ", unique(subdata$conc),
                   ", Nsurv increases at some time points.",
                   sep = "")
      errors <- errorTableAdd(errors, "NsurvIncrease", msg)
    }
    errors
  }
  res <- by(data, list(data$replicate, data$conc), consistency)
  consistency.errors <- do.call("errorTableAppend", res)
  errors <- errorTableAppend(errors, consistency.errors)

  if (diagnosis.plot && "NsurvIncrease" %in% errors$id) {
    survDataPlotFull(data, ylab = "Number of surviving individuals")
  }
  return(errors)
}
