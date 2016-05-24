xslegend <-
function (y, name = NULL)
{
        if (is.factor(y)){
		    par(mar = c(0, 0, 0, 0))
            legend("left", legend = levels(y), 
                fill = factor2color(as.factor(levels(y))),
                title = if (!is.null(name)) name else "", bg = "white")				
		} else {
			par(mar = c(8, 2.2, 8, 2.2))
            fullrange <- abs(diff(range(y)))
			yrange <- seq(min(y, na.rm = TRUE) - 0.15 * fullrange, 
                max(y, na.rm = TRUE) + 0.15 * fullrange, length.out = 80L)
			spacing <- abs(diff(unique(yrange[1:2]))) / 2
            plot(0, 0, xaxt = "n", main = if (!is.null(name)) name else "", 
                ylab = "", col = NULL, pch = 16, xlab = "", bty = "n", 
                xlim = c(-0.5, 0.5), ylim = range(y))
            rect(xleft = -0.7, xright = 0.7, ybottom = yrange - spacing, 
			     ytop = yrange + spacing, col = cont2color(yrange, 
                 range(y)), border = NA)
		    box()
		}
}