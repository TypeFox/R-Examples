quantplot <- 
function(observed, energy = NULL, expected = NULL, bars, 
    barseg = round(length(observed) / bars), 
    main = NULL, xlab = NULL, ylab = "note", xlim = NULL, ylim = NULL, 
    observedcol = "red", expectedcol = "grey", gridcol = "grey",
    lwd = 2, las = 1, cex.axis = 0.9, mar = c(5, 4, 4, 4) + 0.1,
    notenames = NULL, silence = "silence", plotenergy = TRUE, ...,
    axispar = list(ax1 = list(side=1), ax2 = list(side=2), ax4 = list(side=4)),
    boxpar = list(), 
    energylabel = list(text="energy", side=4, line=2.5, at=rg.s-0.25, las=3),
    energypar = list(pch = 20), 
    expectedpar = list(),
    gridpar = list(gridbar = list(col = 1), gridinner = list(col=gridcol)),
    observedpar = list(col = observedcol, pch = 15))
{
    nam_ob <- names(observed)
    if(length(nam_ob) && is.list(observed) && all(c("energy", "notes") %in% nam_ob)){
        if(is.null(energy)) energy <- observed$energy
        observed <- observed$notes
    }
    if(is.null(energy) || all(is.na(energy))) 
        plotenergy <- FALSE
        
    opar <- par(las = las, cex.axis = cex.axis, mar = mar)
    on.exit(par(opar))
  
    if(is.null(xlab)) xlab <- "bar"
    rg <- range(observed, expected, na.rm = TRUE)
    if(is.null(notenames)) notenames <- notenames(rg[1]:rg[2])
    y.ticks <- c(silence, notenames)
    rg.s <- rg[1] - 2
    observed[is.na(observed)] <- rg.s
    x <- seq(along = observed) / barseg
    if(is.null(xlim)) xlim <- c(0, bars)
    if(is.null(ylim)) ylim <- c(rg.s - 2, rg[2] + 0.5)
    
    ## setup pf the plot itself:    
    plot(x, observed, xaxt = "n", yaxt = "n", 
        main = main, xlab = xlab, ylab = ylab, 
        type = "n", lwd = lwd, xaxs = "i", yaxs = "i",
        xlim = xlim, ylim = ylim, ...)
    do.call("axis", c(list(at = 1:bars), axispar[["ax1"]]))
    
    if(!is.null(expected)) 
        do.call("rect", 
            c(list(xleft = c(0, x[-length(x)]), ybottom = expected-.5, 
                   xright = x,                  ytop =    expected+.5),
                   border = expectedcol, col = expectedcol,
              expectedpar))        
    do.call("axis", 
        c(list(at = (rg.s:rg[2])[-2], label = as.character(y.ticks)), 
          axispar[["ax2"]]))
    do.call("points", c(list(x = x - 1 / (2 * barseg), y = observed), observedpar))    

    ## grid:
    do.call("abline", c(list(
        v = 1:(bars*barseg) / barseg, h = rg[1]:rg[2] - 0.5), 
        gridpar[["gridinner"]])) 
    do.call("abline", c(list(v = 1:bars, h = rg.s + 1.5, gridpar[["gridbar"]])))

    if(plotenergy){
        energy1 <- rg.s - 2 + 
            (3.5 * (energy - min(energy)) / diff(range(energy)))
        do.call("points", c(list(
            x = x - 1 / (2 * barseg), 
            y = energy1), energypar))
        do.call("mtext", energylabel)
        do.call("axis", 
            c(list(at = rg.s + c(-2, 1.5), 
                label = round(range(energy), 1)), 
            axispar[["ax4"]]))
    }
    # clean up with a final box:
    do.call("box", boxpar)
}
