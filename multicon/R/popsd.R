popsd <-
function(x, nomiss=.8) {
    miss <- sum(is.na(x))
    comp <- (length(x) - miss) / length(x)
    out <- ifelse(comp >= nomiss, sqrt((var(x, na.rm=T)*(length(x)-miss-1)) / (length(x)-miss)), NA)
    return(out)
}
