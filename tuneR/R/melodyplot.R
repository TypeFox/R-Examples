melodyplot <- 
function(object, observed, expected = NULL, bars = NULL, main = NULL, 
    xlab = NULL, ylab = "note", xlim = NULL, ylim = NULL, 
    observedtype="l", observedcol = "red", expectedcol = "grey", gridcol = "grey",
    lwd = 2, las = 1, cex.axis = 0.9, mar = c(5, 4, 4, 4) + 0.1,
    notenames = NULL, thin = 1, silence = "silence", plotenergy = TRUE, ...,
    axispar = list(ax1 = list(side=1), ax2 = list(side=2), ax4 = list(side=4)),
    boxpar = list(), 
    energylabel = list(text="energy", side=4, line=2.5, at=rg.s-0.25, las=3),
    energypar = list(), 
    expectedpar = list(),
    gridpar = list(col=gridcol), 
    observedpar = list(col=observedcol, type=observedtype, lwd=2, pch=15))
{
    observed <- as.matrix(observed)
    if(!is.null(expected)){
        expected <- as.matrix(expected)
        if(length(dim(observed)) && any(dim(observed) != dim(expected)))
            stop("Dimensions of 'observed' and 'expected' must be equal")
    }
    opar <- par(las = las, cex.axis = cex.axis, mar = mar)
    on.exit(par(opar))

    if(is.null(bars)){
        starts <- object@starts
        bars <- starts[length(starts)] + object@width
        bars <- bars / object@samp.rate 
        if(is.null(xlab)) xlab <- "time"
    }
    else if(is.null(xlab)) xlab <- "bar"
    
    rg <- range(observed, expected, na.rm = TRUE)
    if(!is.null(ylim)){
        if(any(ylim %% 1)) stop("ylim must be a vector of two integers indicating note heights")
        rg <- ylim
    }
    if(is.null(notenames)) notenames <- notenames(rg[1]:rg[2])
    y.ticks <- c(silence, notenames)
    rg.s <- rg[1] - 2
    ylim <- c(rg.s - 2 , rg[2] + 0.5)
    if(is.null(xlim)) xlim <- c(0, if(bars < 1) 1 else bars)

    observed[is.na(observed)] <- NaN
    observed[observed < rg[1] & observed != rg.s] <- NA
    observed[is.nan(observed)] <- rg.s
    if(!is.null(expected)) 
        expected[expected < rg[1] & expected != rg.s] <- NA
    
    x <- bars * seq(0, 1, length = nrow(observed))

    ## setup pf the plot itself:
    plot(x, observed[,1] - .05, xaxt = "n", yaxt = "n", 
        main = main, xlab = xlab, ylab = ylab, 
        type = "n", lwd = lwd, xaxs = "i", yaxs = "i",
        xlim = xlim, ylim = ylim, ...)

    if(!is.null(expected)) 
        do.call("rect", 
            c(list(xleft = c(0, x[-length(x)]), ybottom = expected-.5, 
                   xright = x,                  ytop =    expected+.5),
                   border = expectedcol, col = expectedcol,
              expectedpar))

    ## grid:
    do.call("abline", c(list(v = 1:(bars), h = rg[1]:rg[2] - 0.5), gridpar)) 

    do.call("points", 
        c(list(x=rep(x, dim(observed)[2]), y=observed + .05), observedpar))

    do.call("axis", c(list(at = 1:bars), axispar[["ax1"]]))
    at <- c(rg.s, seq(rg.s+2, rg[2], thin))
    label <- as.character(y.ticks)[c(1, seq(2, length(y.ticks), thin))]
    do.call("axis", 
        c(list(at = at, label = label), 
          axispar[["ax2"]]))

    if(plotenergy){
        ## separates two parts of the plot:
        do.call("abline", c(list(h = rg.s + 1.5), boxpar)) 
        
        ## calculate and plot energy:
        energy <- object@energy
        energy <- rg.s - 2 + 
            (3.5 * (energy - min(energy)) / diff(range(energy)))
        do.call("lines", c(list(x=x, y=energy), energypar))
        do.call("mtext", energylabel)
        do.call("axis", 
            c(list(at = rg.s + c(-2, 1.5), 
                label = round(range(object@energy), 1)), 
            axispar[["ax4"]]))
    }
        
    # clean up with a final box:
    do.call("box", boxpar)
}
