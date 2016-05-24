if (getRversion() >= "2.15.1") globalVariables(c(".valid"))

#' @encoding UTF-8
#' @title Method to Produce Descriptive Statistics Summary
#' @param ... Parameters which are typically ignored
#' @return A data.frame of descriptive statistics
#' @export
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
detail <- function(data...)
  UseMethod("detail")


print.detail<-function(data,...){
  cat("detail: ")
  print(data$call)
}
NULL

#' @encoding UTF-8
#' @title Default Summary Statistics Function
#'
#' @description This function provides up to 14 statistics for an entire data object: number of cases, mean, standard deviation, variance, standard error, median, mad (median absolute deviation), trimmed and winsorized means, range, minimum, maximum, skewness, and kurtosis. Statistics for a factor variable might be computed based on its `levels`, and is shown accompained whit ans \verb{"*"}.
#'
#' @details Trimming is not winsorizing. The winsorization process is more complex than simply excluding data. For example, while in a trimmed estimator the extreme values are discarded, in a winsorized estimator, they are rather replaced by certain percentiles.
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#' @param .data a data object (vector or data.frame).
#' @param by a factor variable
#' @param basic indicates if only a short version of the descriptive table should be returned, the default is \code{basic=TRUE}.
#' @param na.rm a logical value for \code{na.rm}, default is \code{na.rm=TRUE}.
#' @param trim is the proportion of the data to be replaced for estimating the average
#' @param type a numeric value (fraction) to be trimmed. The value in trim will be discarded from the top and bottom of data. See in details below
#' @param k a numeric value for observations in the data set to be discarded while computing the winsorized mean. See details below
#'
#' @return A data frame containing the require computations
#'
#' @keywords Descriptive
#' @keywords Tables
#'
#' @references Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S.}. Springer.
#'
#' @examples
#' #load some data
#' data(marriage)
#'
#' # To apply the function
#' detail(marriage, trim = 0.5, k = 3)
#'
#' @export
detail <-
  function(.data, by = NULL, basic = FALSE, na.rm = TRUE, trim = 0.2, type = 2, k = 1)
  {
	  data_name <- deparse(substitute(.data))
	   s1 <- Sys.time()
    cl <- match.call()
    .valid <- function(.data) {
      sum(!is.na(.data))
    }
    if (!na.rm)
      .data <- stats::na.omit(.data)

    if (is.null(dim(.data)[2])) {
      len <- 1
      stats = matrix(rep(NA, 13), ncol = 13)
      stats[1, 1] <- .valid(.data)
      stats[1, 2] <- mean(.data, na.rm = na.rm)
      stats[1, 3] <- stats::sd(.data, na.rm = na.rm)
      stats[1, 4] <- stats::var(.data, na.rm = na.rm)
      stats[1, 5] <- se(.data, na.rm = na.rm)
      stats[1, 6] <- stats::median(.data, na.rm = na.rm)
      stats[1, 7] <- stats::mad(.data, na.rm = na.rm)
      stats[1, 8] <- mean(.data, na.rm = na.rm, trim = trim)
      stats[1, 9] <- winsorize(.data, k=k, na.rm = na.rm)
      stats[1, 10] <- min(.data, na.rm = na.rm)
      stats[1, 11] <- max(.data, na.rm = na.rm)
      stats[1, 12] <- skewness(.data, na.rm = na.rm, type = type)
      stats[1, 13] <- kurtosis(.data, na.rm = na.rm, type = type)

      vars <- 1
    }
    else {
      stats = matrix(rep(NA, ncol(.data) * 13), ncol = 13)
      rownames(stats) <- colnames(.data)
      stats[, 1] <- apply(.data, 2, .valid)
      vars <- c(1:ncol(.data))
      for (i in seq_along(.data)) {
        if (is.factor(.data[[i]]) || is.logical(.data[[i]])) {
          .data[[i]] <- as.numeric(.data[[i]])
          rownames(stats)[i] <- paste(rownames(stats)[i],
                                      "*", sep = " ")
        }
        .data[sapply(.data, is.character)] <- lapply(.data[sapply(.data, is.character)],
                                             as.factor)
        if (!is.numeric(unclass(.data[[i]])))
          stop("'detail' still doesn't know how to deal with  'Date' or 'character' vectors")
      }
      stats[, 2] <- apply(.data, 2, mean, na.rm = na.rm)
      if (!basic) {
        stats[, 12] <- sapply(.data, FUN=skewness, na.rm = na.rm, type = type)
        stats[, 13] <- sapply(.data, FUN=kurtosis, na.rm = na.rm, type = type)
      }
      stats[, 3] <- sapply(.data, FUN=stats::sd, na.rm = na.rm)
      stats[, 4] <- sapply(.data, FUN=stats::var, na.rm = na.rm)
      stats[, 5] <- sapply(.data, FUN=se, na.rm = na.rm)
      stats[, 6] <- sapply(.data, FUN=stats::median, na.rm = na.rm)
      stats[, 7] <- sapply(.data, FUN=stats::mad, na.rm = na.rm)
      stats[, 8] <- sapply(.data, FUN=mean, na.rm = na.rm, trim = trim)
      stats[, 9] <- sapply(.data, FUN=winsorize, k=k, na.rm = na.rm)
      stats[, 10] <- sapply(.data, FUN=min, na.rm = na.rm)
      stats[, 11] <- sapply(.data, FUN=max, na.rm = na.rm)
    }
      if (!basic) {
        temp <- data.frame(Obs = stats[, 1],
                           Mean = stats[, 2],
                           SD = stats[, 3],
                           Var = stats[, 4],
                           SE = stats[, 5],
                           Median = stats[, 6],
                           Mad = stats[, 7],
                           Trimmed = stats[, 8],
                           winsorize = stats[, 9],
                           Range = stats[, 11] - stats[, 10],
                           Min = stats[, 10],
                           Max = stats[, 11],
                           Skew = stats[, 12],
                           Kurt = stats[, 13])
      }
      else {
        temp <- data.frame(Obs = stats[, 1],
                           Mean = stats[, 2],
                           SD = stats[, 3],
                           Min = stats[, 10],
                           Max = stats[, 11] )
      }

    output <- format(round(data.frame(temp), 1), nsmall = 0)
    class(output) <- c("SciencesPo", "describe", "data.frame")
    return(output)

s2 <- Sys.time()
timediff <- c( s2 - s1 )
cat("\n")
cat("Date of Analysis: ", format(Sys.time(), "%a %b %d %Y"), "\n", "Computation time: ",timediff, sep="", "\n")
cat("\n")
cat("Descriptive Statistics for ", length(vars), " variables", " from ", data_name, sep = "","\n")
cat("-------------------------------------------------------------\n")
cat("\n")
  }
