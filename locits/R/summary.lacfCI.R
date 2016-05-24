summary.lacfCI <-
function (object, ...) 
{
	CI <- object
	object <- object$the.lacf
	nlags <- dim(object$lacf)
	ntime <- nlags[1]
	nlags <- nlags[2]
    cat("Name of originating time series: ", object$name, "\n")
    cat("Date produced: ", object$date, "\n")
    cat("Number of times: ", ntime, "\n")
    cat("Number of lags: ", nlags, "\n")
    cat("Number of CI lags: ", length(CI$lag), "\n")
}
