plotbeta = function(betaestlist, betastderrlist, argvals=NULL,
                    xlab="", ...)
{
#  PLOTBETA plots a functional parameter along with confidence
#  limits
#  Arguments
#  BETAESTLIST    ... A list object containing one or more functional
#                     parameter objects or functional data objects.
#  BETASTDERRLIST ... A list object containing functional data objects
#                     for the standard error of the objects in
#                     BETAESTLIST.

#  Last modified 12 December 2008

#  check BETAESTLIST

if (inherits(betaestlist, "fdPar") || inherits(betaestlist, "fd")) {
    betaestlist = list(betaestlist)
}

if (!inherits(betaestlist, "list")) {
    stop("BETAESTLIST is not a list, fd, or fdpar object.")
}

#  check BETASTDERRLIST

if (inherits(betastderrlist, "fd")) {
    betastderrlist = list(betastderrlist)
}

if (!inherits(betastderrlist, "list")) {
    stop("BETASTDERRLIST is not a list, or fd object.")
}

#  get range

if (is.fdPar(betaestlist[[1]])) {
    rangeval = betaestlist[[1]]$fd$basis$rangeval
} else {
    if (is.fd(betaestlist[[1]])) {
        rangeval = betaestlist[[1]]$basis$rangeval
    } else {
        stop(paste("A list does not contain either a functional parameter ",
           "or a functional data object."))
    }
}

if (is.null(argvals)) {
    argvals = seq(rangeval[1],rangeval[2],len=51)
}
n = length(argvals)
p = length(betaestlist)

par(ask=T)
for (j in 1:p) {
    if (is.fdPar(betaestlist[[j]])) {
        betavec = eval.fd(argvals, betaestlist[[j]]$fd)
    } else {
        if (is.fd(betaestlist[[j]])) {
            betavec = eval.fd(argvals, betaestlist[[j]])
        } else {
            stop(
        "BETAESTLIST does not contain a functional parameter or data object.")
        }
    }
    betastderr = eval.fd(argvals, betastderrlist[[j]])
    betavecp   = betavec + 2*betastderr
    betavecm   = betavec - 2*betastderr
    zeroval  = c(0,0)
    plot(argvals, betavec, type="l", xlab=xlab, ylab="",
         xlim=rangeval, ylim=c(min(betavecm),max(betavecp)), ...)
    lines(rangeval, zeroval,lty=3, col=2)
        #fill([argvals   argvals(n-11)], ...
        #     [betavecp   betavecm(n-11)], "w", "LineStyle", "--")
    lines(argvals, betavec,  col=1, lwd=2)
    lines(argvals, betavecp, col=1, lwd=1)
    lines(argvals, betavecm, col=1, lwd=1)
    lines(rangeval, zeroval, col=1, lty=3)
    for (i in 1:n) {
            lines(c(argvals[i],argvals[i]),c(betavecm[i],betavecp[i]))
    }
    title(paste("Regression function ",j))
}

}
