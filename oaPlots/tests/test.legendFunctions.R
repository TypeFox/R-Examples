library(oaPlots)
library(oaColors)

###################
# layout tests


# 1
layout = matrix(c(1,1,0,2), 2, 2, byrow=TRUE)
type = "layout"
side = "bottom"
heights = c(3, 1)
widths = c(2, 1)
proportion = 0.1
n <- 7


prepLegend(layout = layout, side = side, proportion = proportion, 
		heights = heights, width = widths)
for(i in 1:max(layout))
	plot(1:n, 1:n, col = oaPalette(n), pch = 19, cex = 2.2, xaxt = "n", 
			yaxt = "n", ann = FALSE)
addLegend(legend = names(oaPalette(n)), font = 2, 
		pch = 19, pt.cex = 2.2, text.col = oaPalette(n), col = oaPalette(n), 
		text.width = .11, xjust = 0, horiz = TRUE)

# 2
layout = rbind(c(1, 2, 3), c(0, 4, 3), c(0, 4, 5))
type = "layout"
side = "right"
heights = c(1, 1, 1)
widths = c(1, 1, 1)
proportion = 0.15

prepLegend(layout = layout, side = side, proportion = proportion, type = "layout",
		heights = heights, width = widths)
for(i in 1:max(layout))
	plot(1:n, 1:n, col = oaPalette(n), pch = 19, cex = 2.2, xaxt = "n", 
			yaxt = "n", ann = FALSE)
addLegend(legend = names(oaPalette(n)), font = 2, 
		pch = 19, pt.cex = 2.2, text.col = oaPalette(n), col = oaPalette(n), 
		y.intersp = 1.5, xjust = 0, cex = 1.5)


# 3
layout = rbind(c(1, 2), c(0, 2), c(3, 2))
type = "layout"
side = "left"
heights = NULL
widths = NULL
proportion = 0.15

prepLegend(layout = layout, side = side, proportion = proportion, type = "layout",
		heights = heights, width = widths)
for(i in 1:max(layout))
	plot(1:n, 1:n, col = oaPalette(n), pch = 19, cex = 2.2, xaxt = "n", 
			yaxt = "n", ann = FALSE)
addLegend(legend = names(oaPalette(n)), font = 2, 
		pch = 19, pt.cex = 2.2, text.col = oaPalette(n), col = oaPalette(n), 
		y.intersp = 1.5, xjust = 0, cex = 1.5)




###################
# mfrow tests


# 1
layout <- c(1,1);
side <- "right"
proportion <- 0.15
n <- 7


prepLegend(layout = layout, side = side, proportion = proportion)
for(i in 1:(layout[1]*layout[2]))
	plot(1:n, 1:n, col = oaPalette(n), pch = 19, cex = 2.2, xaxt = "n", 
			yaxt = "n", ann = FALSE)
addLegend(x = "left", legend = names(oaPalette(n)), font = 2, 
		pch = 19, pt.cex = 2, text.col = oaPalette(n), col = oaPalette(n), 
		y.intersp = 1.5, xjust = 0)



# 2
layout <- c(2,3);
side <- "left"
proportion <- 0.2

prepLegend(layout = layout, side = side, proportion = proportion)
for(i in 1:(layout[1]*layout[2]))
	plot(1:7, 1:7, col = oaPalette(7), pch = 19, cex = 2.2, xaxt = "n", 
			yaxt = "n", ann = FALSE)
addLegend(legend = names(oaPalette(7)), font = 2, 
		pch = 19, pt.cex = 2, text.col = oaPalette(7), col = oaPalette(7), 
		y.intersp = 1.5, cex = 1.5)


# 3
layout <- c(2,2);
side <- "bottom"
proportion <- 0.05


prepLegend(layout = layout, side = side, proportion = proportion)
for(i in 1:(layout[1]*layout[2]))
	plot(1:n, 1:n, col = oaPalette(n), pch = 19, cex = 2.2, xaxt = "n", 
			yaxt = "n", main = "Lorem Ipsum", xlab = "", ylab = "")
addLegend(legend = names(oaPalette(n)), font = 2, 
		pch = 19, pt.cex = 2, text.col = oaPalette(n), col = oaPalette(n), 
		y.intersp = 1.5, cex = 1.5, horiz = TRUE)


