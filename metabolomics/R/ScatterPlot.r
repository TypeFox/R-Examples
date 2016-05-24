ScatterPlot <- function(inputdata, y.axis, x.axis, 
    ylab=y.axis, xlab=x.axis, ...)
{
    groups <- levels(factor(inputdata[, 1], levels = unique(inputdata[, 1])))
    # Create a list of point symbols to use
    pch_list <- c(                         # square, circle, triangle, diamond
        15, 16, 17, 18,                    # as solid
        0, 1, 2, 5,                        # as open
        7, 13, 9                           # as open with cross (no triangle)
    )

    # Create a list of colours (chosen to be colourblind safe)
    col_list <- ColList(15)
    
    # Determine whether or not to cycle through colours
    col_test <- ceiling(length(groups) / length(col_list))
    pch_test <- ceiling(length(groups) / length(pch_list))
    
    if (pch_test > 1) {
        col_list <- rep(col_list, col_test)
    }
    if (col_test > 1) {
        pch_list <- rep(pch_list, pch_test)
    }
    
    # Find longest label in legend
    leg_width <- 1
    for (ii in 1:length(groups)) {
        if (nchar(groups[ii]) > leg_width) {
            leg_width <- nchar(groups[ii])
        }
    }
    
    lax_width <- 1
    for (ii in 1:length(inputdata[[x.axis]])) {
        if (nchar(inputdata[[x.axis]][ii]) > lax_width) {
            lax_width <- nchar(inputdata[[x.axis]][ii])
        }
    }
    for (ii in 1:length(inputdata[[y.axis]])) {
        if (nchar(inputdata[[y.axis]][ii]) > lax_width) {
            lax_width <- nchar(inputdata[[y.axis]][ii])
        }
    }

    # Save original graphical parameters
                        ## Can also be written as
    #par_defs <- par(
    #    font=par()$font,
    #    xpd=par()$xpd,
    #    mar=par()$mar,
    #    ...
    #)
    par_defs <- par(
        font=1,
        xpd=FALSE,
        mar=c(5.1, 4.1, 4.1, 2.1), ...
    )
    # Change as desired
    l_mar_exp <- 5.5 * strwidth(lax_width, units="i")
    r_mar_exp <- 5.5 * strwidth(leg_width, units="i") + leg_width + 0.25
    # Set new graphical parameters
    par(font=2,
        xpd=TRUE,
        mar=par()$mar + c(0, l_mar_exp, 0, r_mar_exp)
    )
    
    #
    #    Plotting
    #
    plot(range(inputdata[[x.axis]]),
        range(inputdata[[y.axis]]),
        las=1,                             # axis labels always horizontal
        type="n",                          # don't plot anything (yet)
        xlab=xlab,                         # x-axis title
        ylab=ylab
    )
    
    # Plot points
    for (ii in 1:length(groups)) {
        points(inputdata[which(inputdata[, 1] == groups[ii]), x.axis],
            inputdata[which(inputdata[, 1] == groups[ii]), y.axis],
            pch=pch_list[ii],
            col=col_list[ii]
        )
    }
    
    # Add legend
    legend("left",                         # position of legend
        inset=1.01,                        # 1 is RHS of plot area
        groups,                            # labels on legend
        pch=pch_list,                      # points (see above)
        col=col_list,                      # colour (see above)
        cex=1                              # character size
    )
    
    # Reset graphical parameters to defaults
    par(par_defs)
}
