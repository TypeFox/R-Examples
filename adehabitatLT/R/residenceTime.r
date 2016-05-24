residenceTime <- function(lt, radius, maxt, addinfo=FALSE,
                          units = c("seconds", "hours", "days"))
{
    if (!inherits(lt, "ltraj"))
        stop("lt should be of class ltraj")
    if (length(radius)>1)
        stop("Only one radius allowed in this function")
    units <- match.arg(units)
    if (units == "hours")
        maxt <- maxt*3600
    if (units == "days")
        maxt <- maxt*(3600 * 24)

    res <- lapply(1:length(lt), function(i) {
        x <- lt[[i]]
        uu <- x$date
        vv <- .Call("residtime", x, radius, maxt, PACKAGE="adehabitatLT")
        if (all(is.na(vv)))
            warning(paste("Too large radius for burst", burst(lt)[i],"\n",
                          "The residence time is missing for all the relocations of this burst\n"))
        z <- data.frame(uu,vv)
        names(z) <- c("Date", paste("RT", format(radius, scientific=FALSE), sep="."))
        return(z)
    })
    if (addinfo) {
        res <- lapply(res, function(x) {
            x <- data.frame(x[,2])
            names(x) <- paste("RT", format(radius, scientific=FALSE), sep=".")
            return(x)
        })
        if (!is.null(infolocs(lt))) {
            il <- infolocs(lt)
            res <- lapply(1:length(il), function(j) {
                cbind(il[[j]], res[[j]])
            })
        }
        infolocs(lt) <- res
        res <- lt
    } else {
        names(res) <- burst(lt)
        class(res) <- "resiti"
        attr(res, "radius") <- radius
        attr(res, "maxt") <- maxt
    }
    return(res)
}


print.resiti <- function(x, ...)
{
    cat("*****************************\n")
    cat("* Object of class resiti\n")
    cat("* (residence time method)\n\n")
    cat("Radius =", attr(x, "radius"), "\n")
    cat("Maximum time allowed before coming back in the circle =", attr(x, "maxt"), "seconds\n")
    cat("This object is a list of data.frames.\nThe following bursts are available:\n")
    cat(paste("$", names(x), "\n", sep=""))
}

plot.resiti <- function(x, addpoints=FALSE, addlines=TRUE, ...)
{
    par(mfrow=n2mfrow(length(x)))
    tmp <- lapply(1:length(x), function(i) {
        y <- x[[i]]
        plot(y[,1], y[,2], ty="n", xlab="Date",
             ylab=paste("Residence Time in a circle of", attr(x, "radius")),
             main=names(x)[i])
        if (addpoints) {
            points(y[,1], y[,2], pch=16, ...)
        }
        if (addlines)
            lines(y[,1], y[,2])
    })
}
