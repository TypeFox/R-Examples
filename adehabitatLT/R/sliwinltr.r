
sliwinltr <- function(ltraj, fun, step, type=c("locs","time"),
                      units=c("sec", "min", "hour", "day"),
                      plotit = TRUE, ...)
{
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    if ((!attr(ltraj,"typeII"))&(type == "time"))
        stop("time not recorded for ltraj - use type=\"locs\"")

    type <- match.arg(type)
    units <- match.arg(units)


    if (type=="locs") {
        res <- lapply(ltraj, function(x) {
            if (!is.null(attr(x, "infolocs"))) {
                x <- cbind(attr(x, "infolocs"), x)
            }
            uu <- apply(embed(1:nrow(x),step),1, function(y) {
                return(fun(x[y,]))
            })
            return(uu)
        })
        if (plotit) {
            opar <- par(mfrow=n2mfrow(length(ltraj)))
            on.exit(par(opar))
            lapply(1:length(ltraj), function(i) {
                plot(res[[i]], xlab="date", ylab="smoothed values",
                     ty="l", main=names(ltraj)[i])
                points(res[[i]], pch=16)
            })
        }
    } else {
        if (units=="min")
            step <- step*60
        if (units=="hour")
            step <- step*60*60
        if (units=="day")
            step <- step*60*60*24

        res <- lapply(ltraj, function(x) {
            if (!is.null(attr(x, "infolocs"))) {
                x <- cbind(attr(x, "infolocs"), x)
            }
            dam1 <- x$date - step
            dap1 <- x$date + step
            da <- x$date
            re <- as.data.frame(do.call("rbind", lapply(1:length(da), function(i) {
                xt <- x[da>=dam1[i]&da<dap1[i],]
                re <- fun(xt,...)
                return(c(da[i],re))})))
            class(re[,1]) <- c("POSIXct","POSIXt")
            attr(re[,1], "tzone") <- attr(ltraj[[1]]$date, "tzone")
            names(re) <- c("date","y")
            return(re)
        })

        names(res) <- names(ltraj)
        if (plotit) {
            opar <- par(mfrow=n2mfrow(length(ltraj)))
            on.exit(par(opar))
            lapply(1:length(ltraj), function(i) {
                plot(res[[i]][,1], res[[i]][,2],
                     xlab="date", ylab="smoothed values",
                     ty="l", main=names(ltraj)[i])
                points(res[[i]][,1], res[[i]][,2], pch=16)
            })
        }
    }
    invisible(res)
}






