myPoisson <- function(prob) {
    require(sampling)
    which(as.logical(UPpoisson(prob)))
}