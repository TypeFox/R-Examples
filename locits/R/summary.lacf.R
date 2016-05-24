summary.lacf <-
function (object, ...) 
{
    nlags <- dim(object$lacf)
    ntime <- nlags[1]
    nlags <- nlags[2]
    cat("Name of originating time series: ", object$name, "\n")
    cat("Date produced: ", object$date, "\n")
    cat("Number of times: ", ntime, "\n")
    cat("Number of lags: ", nlags, "\n")
}
