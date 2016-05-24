# plotld.R

plotld <- function(x,  # x is passed directly to princomp
    npc=3,
    main="Loadings",
    lty=1,
    lwd=4 / 1:npc,
    col=gray(0:(npc-1) / npc),
    ylim=range(loadings),
    abs.=FALSE,
    cex=.8,
    ylab=if(abs.) "abs(loading)" else "loading",
    legend.x=NULL,
    legend.y=NULL)
{
    loadings <- princomp(x, scores=FALSE)$loadings
    if(anyNA(loadings))
        stop("NA in loadings")
    if(abs.)
        loadings <- abs(loadings)

    if(!is.numeric(npc) || npc < 1)
        stop("invalid npc")
    npc <- min(npc, ncol(loadings))
    col <- rep(col, npc) # ensure we have enough colors and line types
    lty <- rep(lty, npc)
    lwd <- rep(lwd, npc)

    old.par <- par(no.readonly=TRUE)
    on.exit(par(old.par))
    mar <- par("mar")
    mar[1] = 1.2 * mar[1] # increase bottom margin for space for var names
    mar[3] = .5 * mar[3]  # decrease top margin, we don't need it so big
    par(mar=mar)

    plot(loadings[,1], type="l", ylim=ylim, main=main,
         xaxt="n", xlab="", ylab=ylab, col=col[1], lwd=lwd[1])
    abline(h=0, lty=2) # axis line
    # variable names along the x axis
    mtext(rownames(loadings), side=1, at=1:ncol(loadings), cex=cex, las=2, line=1)
    if(npc > 1) {
        for (i in 2:npc)
            lines(loadings[,i], col=col[i], lwd=lwd[i])
        # figure out legend.x and legend.y
        xjust=0
        if(is.null(legend.x)) {
            legend.x <- ncol(loadings)
            xjust <- 1
        }
        if(is.null(legend.y)) {
            # take a stab at guessing if legend should be on top or bottom
            legend.y <- ylim[2]     # legend on top
            yjust <- 1
            nr <- nrow(loadings)
            min <- min(loadings[(nr-1):nr,1:npc])
            max <- max(loadings[(nr-1):nr,1:npc])
            if(abs(ylim[1] - min) > ylim[2] - max) {
                legend.y <- ylim[1] # legend on bottom
                yjust <- 0
            }
        }
        legend(legend.x, legend.y, legend=paste("pc", 1:npc, sep=""),
               col=col, lwd=lwd, bg="white", xjust=xjust, yjust=yjust, cex=cex)
    }
    invisible(loadings)
}
