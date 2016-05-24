"peaks" <- 
function(x, span = 3, strict = TRUE)
{
    span <- as.integer(span)
    if(span %% 2 != 1) {
        span <- span + 1
        cat("span increased to next odd value: ", span, "\n")
    }
    halfspan <- (span - 1)/2
    a.x <- attributes(x)
    dfp <- inherits(x, "data.frame")
    x <- as.matrix(x)
    dx <- dim(x)
    if(length(wna <- which.na(x))) {
        x.ok <- rep(TRUE, length(x))
        x.ok[wna] <- FALSE
        m <- min(x[x.ok])
        x[!x.ok] <- m - 1
    }
    z <- .Fortran("splus2rpeaks",
        as.double(x),
        as.integer(halfspan),
        as.logical(strict),
        as.integer(dx[1]),
        as.integer(dx[2]),
        double(length(x)))[[6]]
    mode(z) <- "logical"
    dim(z) <- dx
    if(dfp) {
        z <- as.data.frame(z)
    }
    attributes(z) <- a.x
    z
}
