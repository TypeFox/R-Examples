impute.mean <- function(x) {
    calc.mean <- mean(x, na.rm = TRUE)
    x[is.na(x)] <- calc.mean
    return(x)
}