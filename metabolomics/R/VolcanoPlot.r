VolcanoPlot <- function(folds, pvals, cexcutoff=0.7, cexlab=0.5, plimit=0.05,
    fclimit=2, xlab='log2 Fold Change', ylab='-log10 t-Test P-value', 
    main="Volcano Plot", ...)
{
    # Get range for x-vals
    x_min <- (-1.5)
    x_max <- 1.5
    if (min(range(folds, finite = TRUE)) <= x_min) {
        x_min <- min(range(folds, finite = TRUE))
    }
    if (max(range(folds, finite = TRUE)) >= x_max) {
        x_max <- max(range(folds, finite = TRUE))
    }
    x_range <- c(x_min, x_max)
    
    # Get range for y-vals
    y_min <- 0
    y_max <- 2
    if (min(range(-log10(pvals), finite = TRUE)) <= y_min) {
        y_min <- min(range(-log10(pvals), finite = TRUE))
    }
    if (max(range(-log10(pvals), finite = TRUE)) >= y_max) {
        y_max <- max(range(-log10(pvals), finite = TRUE))
    }
    y_range <- c(y_min, y_max)
    
    # Draw the basic (empty) plot
    plot(
        x_range,                           # x-dim
        y_range,                           # y-dim
        type='n',                          # empty plot
        xlab=xlab,                         # x-axis title
        ylab=ylab,                         # y-axis title
        main=main,                         # plot title
        ...
    )
    
    # Annotate plot region
    abline(h=-log10(plimit),               # horizontal line at P=plimit
        col='green',                       # line colour
        lty='44'                           # Dot-dash lengths
    )
    mtext(paste("pval =", plimit),         # Label abline
        side=2,                            # on the left plot edge
        at=-log10(plimit),                 # at P=plimit
        cex=cexcutoff,                     # slightly smaller
        las=1                              # perpendicular to axis
    )
    abline(                                # vertical lines at ±2-fold
        v=c(-log2(fclimit), log2(fclimit)),
        col='violet',
        lty='1343'
    )
    mtext(                                 # Label vertical ablines
        c(paste("-", fclimit, "fold"), paste("+", fclimit, "fold")),
        side=3,                            # on top of graph
        at=c(log2(1 / fclimit), log2(fclimit)),
        cex=cexcutoff,
        las=1
    )
    
    # Plot coloured points based on their values
    for (ii in 1:length(pvals)) {
        # If it's below plimit, we're not overly interested: purple.
        if (-log10(pvals[ii]) > (-log10(plimit))) {
            # Otherwise, more checks;
            # if it's greater than 2-fold decrease: blue
            if (folds[ii] > (-1)) {
                # If it's significant but didn't change much: orange
                if (folds[ii] < 1) {
                    points(folds[ii],
                        -log10(pvals[ii]),
                        col='orange',
                        pch=20
                    )
                    text(folds[ii],        # x-coord
                        -log10(pvals[ii]), # y-coord
                        labels=names(folds)[ii],
                        pos=if(-log10(pvals[ii]) < 0.95 * max(y_range)) {
                            if(folds[ii] < 0.75 * max(x_range)) {
                                4          # right if it's neither
                            } else {
                                2          # left if > 0.75 max(x_range)
                            }
                        } else {
                            1              # bottom if > 0.95 max(y_range)
                        },
                        cex=cexlab         # Size of text
                    )
                # Otherwise, greater than 2-fold increase: red
                } else {
                    points(folds[ii],
                        -log10(pvals[ii]),
                        col='red',
                        pch=20
                    )
                    text(folds[ii],        # x-coord
                        -log10(pvals[ii]), # y-coord
                        labels=names(folds)[ii],
                        # If the point is at the top of the
                        # graph, label goes underneath. If it's
                        # at the far right, put the label on
                        # the left of the point.
                        pos=if(-log10(pvals[ii]) < 0.95 * max(y_range)) {
                            if(folds[ii] < 0.75 * max(x_range)) {
                                4          # right if it's neither
                            } else {
                                2          # left if > 0.75 max(x_range)
                            }
                        } else {
                            1              # bottom if > 0.95 max(y_range)
                        },
                        cex=cexlab         # Size of text
                    )
                }
            # Else it's less than -2-fold decrease: blue
            } else {
                points(folds[ii],
                    -log10(pvals[ii]),
                    col='blue',
                    pch=20
                )
                text(folds[ii],            # x-coord
                    -log10(pvals[ii]),     # y-coord
                    labels=names(folds)[ii],
                    # If the point is at the top of the
                    # graph, label goes underneath. If it's
                    # at the far right, put the label on
                    # the left of the point.
                    pos=if(-log10(pvals[ii]) < 0.95 * max(y_range)) {
                        if(folds[ii] < 0.75 * max(x_range)) {
                            4       # right if it's neither
                        } else {
                            2       # left if > 0.75 max(x_range)
                        }
                    } else {
                        1           # bottom if > 0.95 max(y_range)
                    },
                    cex=cexlab  # Size of text
                )
            }
        # Else P > plimit; not significant: purple
        } else {
            points(folds[ii],
                -log10(pvals[ii]),
                col='purple',
                pch=20
            )
        }
    }
}
