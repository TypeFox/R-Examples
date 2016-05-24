savingby2d <- function (x, y = NULL, method = "default")
{
    if(is.data.frame(x) && ncol(x) > 2L) 
        stop("'x' should have max 2 columns.")
    if (is.null(y) && identical(ncol(x), 2L)){
        y <- x[, 2L]
        x <- x[, 1L]
    }
    x <- if (is.data.frame(x))
        x[, 1]
    else x
    y <- if (is.data.frame(y))
        y[, 1]
    else y
    arefactors <- vapply(list(x, y), is.factor, logical(1L))
    if (all(arefactors)){
        tab <- table(x, y)
        return(sum(tab != 0) / (ncol(tab) * nrow(tab)))
    } else {
        if (any(arefactors)){
            if (is.factor(x)){
                fac <- x
                cont <- y
            } else {
                fac <- y
                cont <- x
            }
        totalarea <- abs(diff(range(cont)))
        weightbyfac <- table(fac) / length(fac)
        lengthbyfac <- vapply(levels(fac),
            function(x) {
            if (length(cont[as.character(fac) == x]) > 1)
                abs(diff(range(cont[as.character(fac) == x])))
            else 0
            }, numeric(1))
        hullarea <- sum(weightbyfac * lengthbyfac)
        return(hullarea / totalarea)
        } else {
            if (identical(method, "default")){
                if (abs(cor(x, y)) > 0.995) 
                    return(0)
                x.scaled <- (x - mean(x)) / sd(x)
                y.scaled <- (y - mean(y)) / sd(y)
                totalarea <- abs(diff(range(x.scaled)) * diff(range(y.scaled)))
                conhull <- chull(x.scaled, y.scaled)
                hullarea <- polygonarea(x.scaled[conhull], y.scaled[conhull])
                return(hullarea / totalarea)
            } else {
                if (method %in% c("Outlying", "Skewed", "Clumpy", "Sparse", 
                    "Striated", "Convex", "Skinny", "Stringy", "Monotonic")){
                    if (requireNamespace("scagnostics", quietly = TRUE)){
                        ratio <- scagnostics::scagnostics.default(x, y)[method]
                        if (method %in% c("Outlying", "Skewed", "Clumpy", 
                            "Sparse", "Striated", "Skinny", "Stringy", 
                            "Monotonic"))
                            ratio <- 1 - ratio
                        return(ratio)
                    } else stop("requires package 'scagnostics'")
                    
                } else {
                    if (identical(method, "DECR")){
                        if (requireNamespace("hdrcde", quietly = TRUE)){
                            o <- hdrcde::hdr.2d(x, y, prob = 0.05)
                            return(sum(o$den$z > o$falpha) / length(o$den$z))
                        } else stop("requires package 'hdrcde'")    
                    } else stop("unknown 'method' specified")
                }
            }
        }
    }
}