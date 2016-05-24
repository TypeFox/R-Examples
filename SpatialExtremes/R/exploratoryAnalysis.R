symbolplot <- function(data, coord, which = "gev", plot.border = NULL, col = c("#FF000080", "#0000FF80"),
                       plot.legend = TRUE, scale = 1){

    if (!(which %in% c("gev", "mean", "median")))
        stop("'which' must be one of 'gev', 'mean' or 'median'")

    add <- FALSE
    if (which == "gev"){
        values <- apply(data, 2, gevmle)
        ref.value <- rowMeans(values)

        if (plot.legend)
            layout(matrix(c(2,1,4,3,6,5), 2), heights = c(0.15, 1))
        else
            par(mfrow = c(1, 3))

        mar <- par("mar")

        for (i in 1:3){
            add <- FALSE
            if (!is.null(plot.border)){
                plot.border(add = add)
                add <- TRUE
            }

            sign.col <- col[2 - (values[i,] >= ref.value[i])]
            radius <- abs(values[i,] - ref.value[i])
            norm.factor <-  scale * min(dist(coord)) / max(radius)
            symbols(coord, circles = radius * norm.factor, add = add, bg = sign.col, inches = FALSE)

            if (plot.legend){
                radius.legend <- pretty(radius, 4)
                radius.legend <- radius.legend[radius.legend > 0]
                n.circles <- length(radius.legend)

                xlim <- par("usr")[1:2]
                ylim <- par("usr")[3:4]
                centers <- cbind(seq(xlim[1], xlim[2], length = n.circles), mean(ylim))

                par(mar = c(0, 5, 0, 5))
                plot(-10^6, xlim = xlim, ylim = ylim, axes = FALSE, bty = "n", xlab = "", ylab = "")
                symbols(centers, circles = radius.legend * norm.factor, inches = FALSE, bty = "n",
                        xaxt = "n", yaxt = "n", xlab = "", ylab = "", add = TRUE)

                for (j in 1:n.circles)
                    text(centers[j,1], mean(ylim) - 1.1 * max(radius.legend) * norm.factor,
                         bquote(phantom(a)%+-%.(radius.legend[j])), pos = 1, offset = 0.75)
                
                par(mar = mar)
            }
        }
    }

    else {
                
        if (which == "mean")
            values <- colMeans(data, na.rm = TRUE)
        else
            values <- apply(data, 2, median, na.rm = TRUE)

        ref.value <- mean(values, na.rm = TRUE)

        if (plot.legend)
            layout(matrix(2:1, 2), heights = c(0.15, 1))
        
        mar <- par("mar")
        if (!is.null(plot.border)){
            plot.border(add = add)
            add <- TRUE
        }

        sign.col <- col[2 - (values >= ref.value)]
        radius <- abs(values - ref.value)
        norm.factor <- scale * min(dist(coord)) / max(radius)

        symbols(coord, circles = radius * norm.factor, add = add, bg = sign.col, inches = FALSE)

        for (j in 1:length(radius))
            text(coord[j,1], coord[j,2], bquote(.(round(radius[j], 1))), col = sign.col)

        if (plot.legend){
            radius.legend <- pretty(radius, 4)
            radius.legend <- radius.legend[radius.legend > 0]
            n.circles <- length(radius.legend)

            xlim <- c(0, diff(par("usr")[1:2]))
            centers <- cbind(seq(xlim[1], xlim[2], length = n.circles), 0)

            par(mar = c(0, 5, 0, 5))
            symbols(centers, circles = radius.legend * norm.factor, inches = FALSE, bty = "n",
                    xaxt = "n", yaxt = "n", xlab = "", ylab = "")

            for (j in 1:n.circles)
                text(centers[j,1], -max(radius.legend) * norm.factor,
                     bquote(phantom(a)%+-%.(radius.legend[j])), pos = 1, offset = 0.75)
                
                par(mar = mar)
            }
    }
}


