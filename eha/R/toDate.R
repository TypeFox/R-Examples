toDate <- function(times){
    if (!is.numeric(times)) stop("Argument must be numeric")
    times * 365.2425 + as.Date("0000-01-01")
}
