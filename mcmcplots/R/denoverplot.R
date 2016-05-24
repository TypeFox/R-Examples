denoverplot <- function(mcmc1, mcmc2, parms=NULL, regex=NULL, random=NULL, ci = NULL, auto.layout=TRUE, legend=TRUE, mar=c(2.0, 2.0, 1.5, 0.25)+0.1, col=mcmcplotsPalette(2), lty=1, plot.title=NULL, main=NULL, greek = FALSE, style=c("gray", "plain"), ...){
    nm1 <- deparse(substitute(mcmc1))
    nm2 <- deparse(substitute(mcmc2))

    gpar.args <- list(...)

    style <- match.arg(style)

    ## Convert to mcmc objects if not already
    if (!(is.mcmc(mcmc1)|is.mcmc.list(mcmc1))) mcmc1 <- as.mcmc(mcmc1)
    if (!(is.mcmc(mcmc2)|is.mcmc.list(mcmc2))) mcmc2 <- as.mcmc(mcmc2)

    lty <- rep(lty, length=2)
    col <- rep(col, length=2)

    ## Find common parameters in both mcmc objects
    parnames <- intersect(varnames(mcmc1), varnames(mcmc2))
    if (length(parnames)==0)
        stop("No matching parameter names in 'mcmc1' and 'mcmc2'.  Sorry, chump.")
    parnames <- parms2plot(parnames, parms, regex, random)
    np <- length(parnames)
    nplots <- np + legend

    if (auto.layout){
        op <- mult.fig(nplots, main=plot.title, mar=mar)$old.par
        on.exit(par(op))
    }
    if (is.null(main)){
        main <- parnames
        names(parnames) <- main
    }
    if (greek) {
      main <- .to.greek(main)
    }
    main <- rep(main, length.out=length(parnames))
    names(main) <- parnames

    for (p in parnames){
        denoverplot1(unlist(mcmc1[, p]), unlist(mcmc2[, p]), col=col, lty=lty, style=style, xlab="", ylab="", ci=ci, main=main[p], gpar=gpar.args)
    }

    if (legend){
        plot(c(0, 1), c(0, 1), type="n", yaxt="n", xaxt="n", xlab="", ylab="", bty="n")
        if (style=="plain"){
            rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4])
        } else {
            rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], border=NA, col=gray(0.85))
        }
        lines(c(0, 0.25), c(0.70, 0.70), col=col[1], lty=lty[1], lwd=3)
        lines(c(0, 0.25), c(0.30, 0.30), col=col[2], lty=lty[2], lwd=3)
        text(0.25, 0.70, labels=nm1, pos=4, cex=1.25)
        text(0.25, 0.30, labels=nm2, pos=4, cex=1.25)
    }
}
