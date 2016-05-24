"runsNAltraj" <- function (x, nrep = 500, plotit = TRUE, ...)
{
    if (!inherits(x,"ltraj"))
        stop("x should be of class 'ltraj'")
    opar <- par(mfrow=n2mfrow(length(x)))
    on.exit(par(opar))
    foo <- function(x) {
        x <- as.numeric(is.na(x[,1]))
        n <- length(x)
        if (sum(x)>2) {
            toto <- .C("runsltr", as.integer(x), as.integer(n), double(nrep+1),
                       as.integer(nrep), PACKAGE = "adehabitatLT")[[3]]
            li <- as.randtest(toto[2:length(toto)], toto[1])
            return(li)
        } else {
            return(NULL)
        }
    }

    uu <- lapply(x, foo)
    names(uu) <- unlist(lapply(x, function(i) attr(i, "burst")))
    if (plotit) {
        lapply(1:length(uu), function(i) if (!is.null(uu[[i]]))
               plot(uu[[i]], main=names(uu)[i], ...))
    }
    invisible(uu)
}

