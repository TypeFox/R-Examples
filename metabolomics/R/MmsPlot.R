MmsPlot <- function(inputdata, variables=c("samples", "metabolites"), 
    refvec=NULL, main="Mean, Median and Standard Deviation", ...)
{
    variables <- match.arg(variables)
    if (variables == "metabolites"){
        groups <- as.character(factor(inputdata[, 1], 
            levels=unique(inputdata[, 1]))
        )
        inputdata <- t(inputdata[, -1])
        colnames(inputdata) <- groups
    }
    # Remove groups for data processing
    mms_data <- inputdata[, -1]
    #if "samples"
    # Take the mean of all the variables for each sample
    Mean <- apply(mms_data, 1, mean, na.rm=TRUE)
    # Take the median of all the variables for each sample
    Median <- apply(mms_data, 1, median, na.rm=TRUE)
    # Take the standard deviation of all the variables for each sample
    StdDev <- apply(mms_data, 1, sd, na.rm=TRUE)
    # Join the data into a matrix for plotting
    mms <- data.frame(cbind(Mean, Median, StdDev))

    # Edit the row names if necessary
    mms <- data.frame(t(editcolnames(t(mms))))  ## ?? why assigned to same?

    #
    #    Generate figure and output to the screen
    #
    # Prepare plot characteristics
    ### These have been selected to be useful for both normal vision and
    ### colourblind vision. The multiple redundancy (colour, line style, shape)
    ### is of particular importance, communicating information without
    ### using the names of colours.
    yrange <- range(mms, na.rm = TRUE)
    xrange <- c(1, length(mms$Mean))
    colour_list <- c("#3366aa",            # blue
        "#009872",                         # green
        "#982187"                          # purple
    )
    pchlist <- c(15, 16, 17)               # square, circle, triangle
    linetype <- c(1:3)                     # solid, dashed, dotted

    #
    #    Generate figure and output to the screen
    #
    #x11()           ## not necessary, removed
    par(las=2,                             # y-axis labels horizontal 
        xpd=TRUE,                          # expand plot region if necessary
        mar=par()$mar+c(3,0,0,4),
        ...
    )
    plot(xrange, yrange,                   # set min/max for x and y axes
        main=main,
        type="n",                          # blank plot
        xlab="",                           # x-axis title
        ylab="",                           # y-axis title
        xaxt="n",                          # suppress the x-axis
        lab=c(10, 10, 7),                  # number of tick marks on
        ...                                # axes (see ?par)
    )                                      
    axis(side=1,                           # on the bottom
        labels=rownames(mms),              # use rownames as labels
        at=1:xrange[2]                     # where to put the labels
    )
    # Draw lines
    for (ii in 1:3) {
        lines(c(1:xrange[2]),              # x-values
            mms[, ii],                     # y-values
            type="o",                      # plot points _o_ver lines
            pch=pchlist[ii],               # points (see above)
            col=colour_list[ii],           # colour (see above)
            lty=linetype[ii],              # line type (see above)
            lwd=2                          # line thickness
        )
    }
    
    if(!is.null(refvec)) {
        lines(c(1:xrange[2]), refvec, col="#ee3333", lwd=2, pch=1, type="o")
        legend(#"topleft",                 # position of legend
            "left",                        # put outside plot region
            inset=1.01, 
            c(colnames(mms), "InStd"),     # labels on legend
            pch=c(pchlist, 1),             # points (see above)
            col=c(colour_list, "#ee3333"), # colour (see above)
            lty=c(linetype, 1),            # line type (see above)
            lwd=2,                         # line thickness (see above)
            cex=1                          # character size
        )

    } else {
        legend("left",                     # position of legend
            inset=1.01,                    # put outside plot region
            colnames(mms),                 # labels on legend
            pch=pchlist,                   # points (see above)
            col=colour_list,               # colour (see above)
            lty=linetype,                  # line type (see above)
            lwd=2,                         # line thickness
            cex=1                          # character size
        )
    }
    
    output<-c()
    output$output<-mms
    
    structure(output, class="results")
}
