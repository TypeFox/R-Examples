#' Lattice theme for Analysis of Biological Data
#' 
#' This theme will help produce plots with color scheme similar to the one used
#' in \emph{Analysis of Biological Data}
#' 
#' 
#' @rdname themes
#' @param bw a logical. Use a grayscale theme instead of color?
#' @param lty line types used for \code{panel.superpose}
#' @return a list that can be used as a \code{lattice} theme.
#' @details
#' \code{theme.abd} and \code{col.abd} are the same function with
#' two names.
#' @author Randall Pruim (\email{rpruim@@calvin.edu})
#' @keywords graphics
#' @aliases col.abd,theme.abd
#' @export
#' @examples
#' 
#' trellis.par.set(theme=col.abd(bw=TRUE))
#' show.settings()
#' trellis.par.set(theme=theme.abd(lty=1))
#' show.settings()
#' 
col.abd <- function (bw = FALSE, lty = 1:7) 
{
    aRed <- colorRampPalette(c("white", "darkred"))(10)[10]
    paleRed <- colorRampPalette(c("white", "darkred"))(10)[4]
    lightBlue <- colorRampPalette(c("white", "steelblue"))(10)[5]
    lightRed <- colorRampPalette(c("white", "darkred"))(10)[6]
    veryLightRed <- colorRampPalette(c("white", "darkred"))(12)[3]
    paleGreen <- colorRampPalette(c("white", "darkGreen"))(10)[8]
    if (bw) {
        return(list(
			background = list(col = "transparent"), 
			axis.line = list(col = "gray30"), 
            axis.text = list(col = "gray30"), 
			plot.polygon = list(col = "gray80"), 
            box.rectangle = list(col = "gray10"), 
			box.umbrella = list(col = "gray10", lty = 1), 
			box.dot = list(col = "gray10"), 
			dot.line = list(col = "gray50"), 
            dot.symbol = list(col = "gray30", pch = 16), 
			plot.line = list(col = "black", lwd = 2), 
			plot.symbol = list(col = "black", fill = "gray80", pch = 16), 
			regions = list(col = gray((1:100)/100)), 
            reference.line = list(col = "gray50"), add.line = list(lty = 1, col = "gray80", lwd = 2), 
			superpose.polygon = list(col = c("gray30", "gray70", "black", "gray50", "gray20", 
				"gray80", "gray60", "gray40"), 
			fill = c("gray80")), 
			superpose.line = list(lty = lty, lwd = 2, 
				col = c("gray30", "gray70", "black", "gray50", "gray20", "gray80", "gray60", "gray40")), 
            superpose.symbol = list(
				pch = c(16, 15, 18, 1, 3, 6, 0, 5), 
				cex = rep(0.7, 7), 
				col = c("gray30", "gray70", "black", "gray50", "gray20", "gray80", "gray60", "gray40")), 
			strip.background = list(alpha = 1, col = c("gray80", "gray65")), 
			strip.shingle = list(alpha = 1, col = c("gray60", "gray30")), 
			par.strip.text = list(cex = 0.5)
			))
    }
    else {
        return(list(
			background = list(col = "transparent"), 
			plot.polygon = list(col = "darkred"), 
            box.rectangle = list(col = "darkred"), 
			box.umbrella = list(col = "darkred"), 
            dot.line = list(col = "#e8e8e8"), 
			dot.symbol = list(col = "darkred", pch = 16), 
			plot.line = list(lwd = 2, col = 'darkred'), 
			plot.symbol = list(col = "darkred", pch = 16), 
			# regions = list(col = heat.colors(100)), 
    		regions = list(col= colorRampPalette(c("red", "darkorange", "orange", "yellow", "white"))(100)),
            reference.line = list(col = "#e8e8e8"), 
			add.line = list(lty = 1, col = "gray20", lwd = 2), 
			superpose.polygon = list(
				col = c("darkred", "orange","lightskyblue3", "darkgreen", "darkorchid4", "pink", "lightgreen")
				), 
			superpose.line = list(lty = lty, lwd = 2, 
				col = c("darkred", "orange", "navy", "darkgreen", "lightskyblue3", "darkorchid4", "pink", "lightgreen")
				), 
            superpose.symbol = list(
				pch = c(16, 15, 18, 1, 3, 6, 0, 5), 
				cex = rep(0.7, 7), 
				col = c("darkred", "orange", "navy", "darkgreen", "lightskyblue3", "darkorchid4", "pink", "lightgreen")
				),
			strip.background = list(alpha = 1, 
				col = c(veryLightRed, "#ffe5cc", "#cce6ff", "#ccffff", "#ffccff", "#ffcccc", "#ffffcc")
				), 
            strip.shingle = list(alpha = 1, 
				col = c(lightRed,"#ff7f00", "#0080ff", "#00ffff", "#ff00ff", "#ff0000", "#ffff00")), 
				par.strip.text = list(cex = 0.5)
			))
    }
}

#' @rdname themes
#' @export
theme.abd <- col.abd
