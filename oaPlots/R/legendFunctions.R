
#' Function for arranging plotting layout to accomodate a legend panel
#' @param layout layout vector or matrix
#' @param type type of layout; either "mfrow" or "layout"
#' @param side side of the plot to place legend on; one of "top", "bottom", "left" or "right"
#' @param proportion proportion of plotting window to allocate to legend
#' @param heights height vector for original layout (before the legend panel is appended)
#' @param widths width vector for original layout (before the legend panel is appended)
#' @return none; layout is passed to current device
#' @author Jason Waddell
#' @export
#' @examples
#'layout <- c(2,3);
#'side <- "left"
#'proportion <- 0.2
#'
#'prepLegend(layout = layout, side = side, proportion = proportion)
#'for(i in 1:(layout[1]*layout[2]))
#'	plot(1:7, 1:7, col = 1:7, pch = 19, cex = 2.2, xaxt = "n", 
#'			yaxt = "n", ann = FALSE)
#'addLegend(legend = paste("Group", 1:7), font = 2, 
#'		pch = 19, pt.cex = 2, text.col = 1:7, col = 1:7, 
#'		y.intersp = 1.5, cex = 1.5)
#'
#'
#'layout = rbind(c(1, 2, 3), c(0, 4, 3), c(0, 4, 5))
#'side = "right"
#'proportion = 0.15
#'
#'prepLegend(layout = layout, side = side, proportion = proportion)
#'for(i in 1:max(layout))
#'	plot(1:7, 1:7, col = 1:7, pch = 19, cex = 2.2, xaxt = "n", 
#'			yaxt = "n", xlab = "", ylab = "", main = paste("Plot", i))
#'addLegend(legend = paste("Group", 1:7), font = 2, 
#'		pch = 19, pt.cex = 2, text.col = 1:7, col = 1:7, 
#'		y.intersp = 1.5, cex = 1.5)
prepLegend <- function(layout = c(1,1), type = if(is.matrix(layout)) "layout" else "mfrow",
		side = "right", proportion = 0.15, heights = NULL, widths = NULL){
	
  
  
	if(proportion <= 0 | proportion >= 1)
		stop("The specified 'proportion' must be between 0 and 1.")
	
	if(!(side == "bottom" | side == "top" | 
				side == "left" | side == "right" ) )
		stop("The specified 'position' must be either 'top', 'bottom', 
						'left' or 'right'. ")
	
	if(type == "mfrow" & length(layout) != 2)
		stop("For 'mfrow' specifications, 'layout' must be a vector of length 2.")
	
	
	
	
	
	if(type == "mfrow"){
		layoutMat <- matrix(1:(layout[1]*layout[2]), nrow = layout[1], ncol = layout[2], 
				byrow = TRUE)
		
		switch(side, 
				bottom = {
					layoutMat <- rbind(layoutMat, c(rep(max(layoutMat)+1, layout[2])))
					layout(layoutMat, heights = c(rep(1, layout[1]), 
									layout[1]/(1/proportion - 1)      ))	
				},
				
				top = {
					layoutMat <- rbind(c(rep(max(layoutMat)+1, layout[2])), layoutMat)
					layout(layoutMat, heights = c(layout[1]/(1/proportion - 1), 
									rep(1, layout[1])     ))	
				},
				
				left = {
					layoutMat <- cbind(c(rep(max(layoutMat)+1, layout[1])), layoutMat)
					layout(layoutMat, widths = c(layout[2]/(1/proportion - 1), 
									rep(1, layout[2]))  )
				}, 
				
				right = {
					layoutMat <- cbind(layoutMat, c(rep(max(layoutMat)+1, layout[1])))
					layout(layoutMat, widths = c(rep(1, layout[2]), 
									layout[2]/(1/proportion - 1))     )
				}
		) # close switch(side)
		
		
	} else { # type == "layout"
		
		nCol <- ncol(layout)
		nRow <- nrow(layout)
		
		if(is.null(widths))
			widths = rep(1, nCol)
		
		if(is.null(heights))
			heights = rep(1, nRow)
		
		switch(side, 
				bottom = {
					layoutMat <- rbind(layout, c(rep(max(layout)+1, nCol)))
					layout(layoutMat, widths = widths,
							heights = c(heights, sum(heights)/(1/proportion - 1) )   )
				},
				
				top = {
					layoutMat <- rbind(c(rep(max(layout)+1, nCol)), layout)
					layout(layoutMat, widths = widths,
							heights = c(sum(heights)/(1/proportion - 1), heights)   )
				},
				
				left = {
					layoutMat <- cbind(c(rep(max(layout)+1, nRow)), layout)
					layout(layoutMat, heights = heights, 
							widths = c(sum(widths)/(1/proportion - 1), widths) )
				}, 
				
				right = {
					layoutMat <- cbind(layout, c(rep(max(layout)+1, nRow)))
					layout(layoutMat, heights = heights, 
							widths = c(widths, sum(widths)/(1/proportion - 1)) )
				}
		) # close switch(side)
		
	} # close type == "layout"
	
}



#' Function for adding a legend to an existing device
#' @param x legend x location 
#' @param y legend y location
#' @param legend vector of legend labels
#' @param font legend text font
#' @param bty A character string which determined the type of box which is drawn about plots. If bty is one of "o" (the default), "l", "7", "c", "u", or "]" the resulting box resembles the corresponding upper case letter. A value of "n" suppresses the box.
#' @param xjust how the legend is to be justified relative to the legend x location. A value of 0 means left justified, 0.5 means centered and 1 means right justified.
#' @param yjust the same as xjust for the legend y location.
#' @param ... additional optional arguments to be passed to legend()
#' @return none; legend is added to the current device
#' @author Jason Waddell
#' @export
#' @examples
#'layout <- c(2,3);
#'side <- "left"
#'proportion <- 0.2
#'
#'prepLegend(layout = layout, side = side, proportion = proportion)
#'for(i in 1:(layout[1]*layout[2]))
#'	plot(1:7, 1:7, col = 1:7, pch = 19, cex = 2.2, xaxt = "n", 
#'			yaxt = "n", ann = FALSE)
#'addLegend(legend = paste("Group", 1:7), font = 2, 
#'		pch = 19, pt.cex = 2, text.col = 1:7, col = 1:7, 
#'		y.intersp = 1.5, cex = 1.5)
#'
#'
#'layout = rbind(c(1, 2, 3), c(0, 4, 3), c(0, 4, 5))
#'side = "right"
#'proportion = 0.15
#'
#'prepLegend(layout = layout, side = side, proportion = proportion)
#'for(i in 1:max(layout))
#'	plot(1:7, 1:7, col = 1:7, pch = 19, cex = 2.2, xaxt = "n", 
#'			yaxt = "n", xlab = "", ylab = "", main = paste("Plot", i))
#'addLegend(legend = paste("Group", 1:7), font = 2, 
#'		pch = 19, pt.cex = 2, text.col = 1:7, col = 1:7, 
#'		y.intersp = 1.5, cex = 1.5)
addLegend <- function(x = "center", y = NULL,  legend, font = NULL, bty = "n", 
		xjust = 0.5, yjust = 0.5,  ...){
	
	op <- par(c("mar", "plt"))
	par(mar = c(0,0,0,0))
	par(plt = c(0,1,0,1))
	plot(x = 1, y = 1, type = "n", xlim = c(0,1), ylim = c(0,1),
			bty = "n", axes = FALSE, ann = FALSE)
	
	if(!is.null(font))
		par(font = font)
	
	legend(x = x, y = y, legend = legend, xjust = xjust, yjust = yjust, bty = bty, 
			...)
	par(op)
  
  layout(1) # restore
	
}




