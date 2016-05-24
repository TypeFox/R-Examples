"na.include" <-
function (x) 
{
    if (inherits(x, "data.frame")) {
        for (i in seq(along = x)) x[[i]] <- na.include(x[[i]])
        x
    }
    else {
        if (!((inherits(x, "factor") || !is.null(attr(x, "levels"))) && 
            any(pos <- is.na(x)))) 
            return(x)
        a <- attributes(x)
        a$levels <- l <- c(a$levels, "NA")
        xx <- as.vector(unclass(x))
        xx[pos] <- length(l)
        attributes(xx) <- a
        xx
    }
}
