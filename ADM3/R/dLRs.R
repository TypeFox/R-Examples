`dLRs` <- function(x) {
return(IQR(diff(na.omit(x))) / (4 * qnorm((1 + 0.5) / 2) / sqrt(2)))
}
