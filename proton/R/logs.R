#' The history of logs into the Proton server
#'
#' The dataset describing the history of logs: who, from where and when logged into the Proton server.
#' The subsequent columns in this dataset describe:
#' \itemize{
#'   \item login. The login of the user which logs into the Proton server. 
#'   \item host. The IP address of the computer, from which the log into the Proton server was detected.
#'   \item date. The date of log into the Proton server. Rows are sorted by this column.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name logs
#' @usage data(logs)
#' @format a data frame with 59366 rows and 3 columns.
NULL
