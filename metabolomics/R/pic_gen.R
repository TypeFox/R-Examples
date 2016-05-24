pic_gen <- function(input_matrix, plot_title, 
    plot_labels=rownames(input_matrix), x_label="", y_label="", cex_val=0.7,
    text_on=TRUE, dev='x11', filename='', cols_used, pch_used=NULL)
{
    if (dev == 'x11') {
        dev.new()
    } else if (dev == 'png') {
        png(filename,
            bg="white",                    # background colour
            res=300,                       # image resolution (dpi)
            units="in",
            width=8.3,
            height=5.8                     # image dimensions (inches)
        )
    } else if (dev == 'jpg') {
        jpeg(filename,
            bg="white",                    # background colour
            res=300,                       # image resolution (dpi)
            units="in",
            width=8.3,
            height=5.8,                    # image dimensions (inches)
            quality=100
        )
    } else if (dev == 'tif') {
        tiff(filename,
            bg="white",                    # background colour
            res=300,                       # image resolution (dpi)
            units="in",
            width=8.3,
            height=5.8,                    # image dimensions (inches)
            compression='none'
        )
    }
    
    par(mgp=c(5, 2, 0),                    # axis margins
                                           # (title, labels, line)
        mar=c(6, 6, 4, 2),                 # plot margins (b, l, t, r)
        las=1                              # horizontal labels
    )

    #pch=... for point shape by group - not yet implemented
    # for (grp in groups) {
    # submat <- input_matrix[which(input_matrix[1,] == grp), ]
    plot(input_matrix[, 1],                # x values (PC1)
        input_matrix[, 2],                 # y values (PC2)
        main=plot_title,                   # title of plot
        xlab=x_label,                      # x-axis title
        ylab=y_label,                      # y-axis title
        cex=cex_val,                       # font size
        col=cols_used,                     # coloured points
        pch=pch_used                       # shape of points
    )
    abline(h=0, col="red", lwd=1)          # line across plot area
    abline(v=0, col="blue", lwd=1)
    
    pos_list <- c()
    if (text_on == TRUE) {
        # Label the points
        # Labels on RHS of points if not near top/right edge
        # (pos: bottom=1, left=2, top=3, right=4)
        for (ii in 1:dim(input_matrix)[1]) {
            if(input_matrix[ii, 2] < 0.95 * max(input_matrix[, 2])) {
                if(input_matrix[ii, 1] < 0.75 * max(input_matrix[, 1])) {
                    pos_list[ii] <- 4      # right if it's neither
                } else {
                    pos_list[ii] <- 2      # left if > 0.75 max(x_range)
                }
            } else {
                pos_list[ii] <- 1          # bottom if > 0.95 max(y_range)
            }
        }
        text(input_matrix[, 1],
            input_matrix[, 2],
            plot_labels,
            cex=cex_val,
            pos=pos_list,
            col=cols_used
        )
    }
    
    # Turn off device unless it's onscreen
    if (dev != 'x11') {
        dev.off()
    }
}
