x <- list(createTestObj(return.types=TRUE), as.character(1:16))
objSizes <- matrix(0, nrow=length(x[[1]]), ncol=length(x[[2]]), dimnames=x)
objLens <- objSizes
for (s in seq(ncol(objSizes))) {
    cat("Doing scale=", s, ": ", sep="")
    for (i in seq(nrow(objSizes))) {
        cat("", i)
        x <- createTestObj(i, s)
        objSizes[i, s] <- object.size(x)
        objLens[i, s] <- length(x)
    }
    cat("\n")
}
objSizes
objLens
