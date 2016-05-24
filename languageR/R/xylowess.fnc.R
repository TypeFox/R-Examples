`xylowess.fnc` <-
function(fmla, 
    data, span=2/3,
	  symbolcolor = "darkgrey", linecolor="blue",
    xlabel="", 
    ylabel="",...) {

    # get sensible default labels from the formula
    y = strsplit(as.character(fmla), "[|]")
    if (xlabel == "") xlabel = y[[3]][1]
    if (ylabel == "") ylabel = y[[2]]

    require("lattice", quietly = TRUE)

    # make the xyplot, this code is a simplified version of one of the
    # examples in the on-line help for xyplot
    xyplot(fmla, 
        data = data,
        xlab = xlabel, ylab = ylabel, ...,
        panel = function(x, y) {
            panel.grid(h = -1, v = -1)  
            # only clearly visible on the screen, not on paper
            panel.xyplot(x, y, col.symbol=symbolcolor)
            panel.loess(x, y, span=span, col.line=linecolor)
        }
    )
}

