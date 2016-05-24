plot.Torgegram <- function(x, sp.relationship = c("fc", "fu"),
                           min.cex = 1.5, max.cex = 6, leg.auto = TRUE, main = "", ylab = "",
                           xlab = "Stream Distance", ... )
{
    par.orig <- par(no.readonly = TRUE)
    if(is.null(sp.relationship) | any(is.na(sp.relationship)) |
       length(sp.relationship) > 2) return("sp.relationship mis-specified")
    if(length(sp.relationship)  == 2 & any(sp.relationship != c("fc", "fu")))
        return("sp.relationship mis-specified")
    if(length(sp.relationship)  == 1)
	if(sp.relationship != "fc" & sp.relationship != "fu")
            return("sp.relationship mis-specified")
    if(class(x) != "Torgegram") return(
            "Not a Torgegram object")
    ev <- x

    if(is.null(ev$call$EmpVarMeth)) ev$call$EmpVarMeth <- "MethMoment"
    if(is.null(ev$call$nlag)) ev$call$nlag <- 6
    if(is.null(ev$call$inc)) ev$call$inc <- 0
    if(is.null(ev$call$nlagcutoff)) ev$call$nlagcutoff <- 6
    if(is.null(as.list(match.call()[-1])$ylab)) {
        if(ev$call$EmpVarMeth == "Covariance") {
            ylab <- "Covariance"} else
        ylab <- "Semivariance"
    }
    if(is.null(as.list(match.call()[-1])$col)) {
        colr = c("blue","green")
    }
    else {
        if(length(sp.relationship) == 2 & length(as.list(
                 match.call()[-1])$col) == 1) {
            colr <- rep(as.character(as.list(match.call()[-1])$col),
                        times = 2)
        } else {
            colr <- as.character(as.list(match.call()[-1])$col)
            if(length(colr)>1) colr <- colr[-1]
        }
    }
    if(is.null(as.list(match.call()[-1])$pch)) {
        plch = c(19,19)
    }
    else {
        if(length(sp.relationship) == 2 & length(as.list(
                 match.call()[-1])$pch) == 1) {
            plch <- rep(as.integer(as.character(as.list(match.call()[-1])$pch)),
                        times = 2)
        } else {
            plch <- as.character(as.list(match.call()[-1])$pch)
            if(length(plch)>1) plch <- plch[-1]
            plch <- as.integer(plch)
        }
    }
    if(length(sp.relationship)==2){
        if(is.null(as.list(match.call()[-1])$main)) {
            main <- paste("Estimation Method:", ev$call$EmpVarMeth)
        } else {
            main <- as.list(match.call()[-1])$main
        }
        plot(c(0, max(ev$distance.connect,ev$distance.unconnect)),
             c(min(0, min(ev$gam.connect,ev$gam.unconnect)),
               max(ev$gam.connect,ev$gam.unconnect)),
             type = "n",
             xlab = xlab,
             ylab = ylab,
             main = main, ...)
        maxnp <- max(ev$np.connect,ev$np.unconnect)
        minnp <- min(ev$np.connect,ev$np.unconnect)
        np.range <- maxnp - minnp
        cex.inc <- max.cex - min.cex
        for(i in 1:length(ev$np.connect)) {
            points(ev$distance.connect[i],
                   ev$gam.connect[i], pch = plch[1], col = colr[1],
                   cex = min.cex + cex.inc*ev$np.connect[i]/maxnp)
        }
        for(i in 1:length(ev$np.unconnect)) {
            points(ev$distance.unconnect[i],
                   ev$gam.unconnect[i], pch = plch[2], col = colr[2],
                   cex = min.cex + cex.inc*ev$np.unconnect[i]/maxnp)
        }
        if(leg.auto)
            legend(x = min(ev$distance.connect,ev$distance.unconnect),
                   y = max(ev$gam.connect,ev$gam.unconnect),
                   legend = c("Flow-connected", "Flow-unconnected"), bty = "n",
                   pch = plch,
                   col = c(colr[1], colr[2]))
    } else if(length(sp.relationship)==1){
        if(sp.relationship == "fc") {
            if(is.null(as.list(match.call()[-1])$main)) {
                main <- paste("Flow-connected,  Estimation Method:",
                              ev$call$EmpVarMeth)
            } else {
                main <- as.list(match.call()[-1])$main
            }
            plot(c(0, max(ev$distance.connect)), c(0, max(ev$gam.connect)),
                 type = "n", xlab = xlab, ylab = ylab,
                 main = main, ...)
            nlag <- length(ev$distance.connect)
            maxnp <- max(ev$np.connect)
            minnp <- min(ev$np.connect)
            np.range <- maxnp - minnp
            cex.inc <- max.cex - min.cex
            for(i in 1:nlag) {
                points(ev$distance.connect[i],
                       ev$gam.connect[i], pch = plch[1],
                       cex = min.cex + cex.inc*ev$np.connect[i]/maxnp,
                       col = colr[1])
            }
        } else if(sp.relationship == "fu") {
            if(is.null(as.list(match.call()[-1])$main)) {
                main <- paste("Flow-unconnected,  Estimation Method:",
                              ev$call$EmpVarMeth)
            } else {
                main <- as.list(match.call()[-1])$main
            }
            plot(c(0, max(ev$distance.unconnect)), c(0, max(ev$gam.unconnect)),
                 type = "n", xlab = xlab, ylab = ylab,
                 main = main, ...)
            nlag <- length(ev$distance.unconnect)
            maxnp <- max(ev$np.unconnect)
            minnp <- min(ev$np.unconnect)
            np.range <- maxnp - minnp
            cex.inc <- max.cex - min.cex
            for(i in 1:nlag) {
                points(ev$distance.unconnect[i],
                       ev$gam.unconnect[i], pch = plch[1],
                       cex = min.cex + cex.inc*ev$np.unconnect[i]/maxnp,
                       col = colr[1])
            }
        } else return("sp.relationship mis-specified")
    }
    par(par.orig)
}

