#' Pulls Data from the Thinknum API
#'
#' For details on the api being implemented here, go to www.thinknum.com/API
#' @param expression Thinknum expression specified as a string.
#' @return return value is a data.frame, the first column is a date and all other columns are numerics.
#' @references This R package uses the Thinknum API. For more information go to http://www.thinknum.com/api. 
#' @author Gregory Ugwi
#' @examples \dontrun{
#' googRev = Thinknum("total_revenue(goog)", type="ts")
#' }
#' @importFrom RCurl getURL
#' @importFrom RJSONIO fromJSON
#' @export
Thinknum <- function(expression) {

    ## Build API URL and add auth_token if available
    string <- paste("http://data.thinknum.com/api/v1/?expression=", URLencode(expression, reserved = TRUE), sep="")
    ## Download and parse data
    response <- getURL(string)
    if (length(grep("403 Forbidden", response)))
        stop("Your usage has been flagged as a violation of Thinknum's terms of service agreement. For help, please contact us at support@thinknum.com")
    json <- try(fromJSON(response, nullValue = as.numeric(NA)), silent = FALSE)
    ## Check if code exists
    if (inherits(json, 'try-error')) {
        error_str <- paste("Something has gone wrong. Please copy output and email to thinknum@thinknum.com:", string, sep="\n")
        stop(error_str)
    }
    if (json["error"] == "Requested expression does not exist.")
        stop("Expression does not exist")
    if (length(json$data) == 0)
        stop("Requested Expression does not exist.")

    ## Shell data from JSON's list
    data        <- as.data.frame(matrix(unlist(json$data), ncol = length(json$column_names), byrow = TRUE),stringsAsFactors=FALSE)
    names(data) <- json$column_names
    data[,1]    <- as.Date(as.POSIXct(data[, 1]/1000, origin='1970-01-01', tz='UTC'))

    ## Transform values to numeric
    if (ncol(data) > 2)
        data[, 2:ncol(data)] <- apply(data[, 2:ncol(data)], 2, as.numeric)
    else
        data[, 2] <- as.numeric(data[, 2])
    return(data)
}
