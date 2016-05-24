list2table3l = function(l1, l2, l3){
	numbers = c(length(intersect(setdiff(l1, l2), setdiff(l1, l3))),
	            length(intersect(setdiff(l2, l1), setdiff(l2, l3))),
	            length(intersect(setdiff(l3, l1), setdiff(l3, l2))),
							length(setdiff(intersect(l1, l2), l3)),
							length(setdiff(intersect(l1, l3), l2)),
							length(setdiff(intersect(l3, l2), l1)),
							length(intersect(intersect(l1, l2), l3)))
	percentage = round(numbers / sum(numbers), 2) * 100
	return(data.frame(numbers = numbers, percentage = percentage))
}

list2table2l = function(l1, l2){
	numbers = c(length(setdiff(l1, l2)), length(setdiff(l2, l1)), length(intersect(l1, l2)))
	percentage = round(numbers / sum(numbers), 2) * 100
	return(data.frame(numbers = numbers, percentage = percentage))
}

scale_radius = function(data, const = 0.25){
	data$radius = sqrt(data$percentage)/10 * const
	return(data)
}

bubbleplot2l = function(lists, percentage, colors, scale){
	x = (c(0.20, 0.80, 0.50) - 0.5) * scale + 0.5
	y = c(0.5, 0.5, 0.5)
	data = list2table2l(lists[[1]], lists[[2]])
	data = data.frame(x, y, data)
	data = scale_radius(data, const = 0.23)
	data$colors = colors
	data = data[order(data$numbers, decreasing = TRUE), ]
	
	grid.lines(x = data$x, y = data$y, gp = gpar(col = "grey60"))
	grid.circle(x = data$x, y = data$y, r = data$radius, gp = gpar(fill = data$colors, col = "grey60"))
	if(percentage){
		grid.text(label = data$percentage, x = data$x, y = data$y)
	}
	else{
		grid.text(label = data$numbers, x = data$x, y = data$y)
	}
	
	grid.text(label = names(lists), x = x[c(1, 2)], y = y[c(1, 2)] + max(data$radius) * 1.05, vjust = 0, hjust = c(0.5, 0.5), gp = gpar(cex = 1.3))
}

bubbleplot3l = function(lists, percentage, colors, scale){
	data = list2table3l(lists[[1]], lists[[2]], lists[[3]])
	x = (c(0.20, 0.50, 0.80, 0.35, 0.50, 0.65, 0.50) - 0.5) * scale + 0.5
	y = (c(0.800, 0.284, 0.800, 0.536, 0.800, 0.536, 0.632) - 0.5) * scale + 0.5
	data = data.frame(x, y, data)
	data = scale_radius(data, const = 0.23)
	data$colors = colors
	
	triangle =  data[c(1, 2, 3, 1, 7, 3, 2, 7), c("x", "y")]
	data = data[order(data$numbers, decreasing = TRUE), ]
	
	grid.lines(x = triangle$x, y = triangle$y, gp = gpar(col = "grey60"))
	grid.circle(x = data$x, y = data$y, r = data$radius, gp = gpar(fill = data$colors, col = "grey60"))
	if(percentage){
		grid.text(label = data$percentage, x = data$x, y = data$y)
	}
	else{
		grid.text(label = data$numbers, x = data$x, y = data$y)
	}
	grid.text(label = names(lists), x = x[c(1:3)], y = y[c(1:3)] - c(-1, 1, -1) * max(data$radius) * 1.05, vjust = c(0, 1, 0), gp = gpar(cex = 1.3))
}

vplayout = function(x, y){
	return(viewport(layout.pos.row = x, layout.pos.col = y))
}

 
#' Simple alternative to Venn diagrams
#' 
#' A simple alternative to the traditional Venn diagram. It depicts each overlap as a 
#' separate bubble with area proportional to the overlap size.
#' 
#' Colors can be specified as vector. For 2 set case a 3 element vector is required with 
#' colors for: Set1, Set2 and Set1 & Set2 correspondingly. For 3 set case a 7 element 
#' vector is required with colors for: Set1, Set2, Set3, Set1 & Set2, Set1 & Set3, Set2 & 
#' Set3 and Set1 & Set2 & Set3 correspondingly.
#'
#' @param sets list of vectors to overlap. If list contains more than 3 elements only the 
#' first 3 are used.
#' @param percentage logical showing if percentages or raw numbers are displayed
#' @param colors vector of colors for the bubbles, see details on specifying that
#' @param fontsize fontsize used for the numbers in the bubbles
#' @param main title of the plot
#' @param scale a scaling factor to adjust the base triangle size when the plot 
#' does not fit the window well. 
#' @param add logical determining if the figure is added to exixting plot or if a new 
#' plot is initialized
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' @examples
#' 
#' bvenn(list(Set1 = sample(letters, 14), Set2 = sample(letters, 9)))
#' bvenn(list(Set1 = sample(letters, 16), Set2 = sample(letters, 12), Set3 = sample(letters, 7)))
#' 
#' # Adding colors
#' bvenn(list(Set1 = sample(letters, 14), Set2 = sample(letters, 9)), colors = c("red", 
#' "green", "yellow"))
#' bvenn(list(Set1 = sample(letters, 16), Set2 = sample(letters, 12), Set3 = 
#' sample(letters, 7)), colors = c("red", "blue", "yellow", "purple", "orange", "green", 
#' "brown"))
#' 
#' # Adjust the triangle size
#' bvenn(list(Set1 = sample(letters, 16), Set2 = sample(letters, 12), Set3 = 
#' sample(letters, 7)), colors = c("red", "blue", "yellow", "purple", "orange", "green", 
#' "brown"), scale = 0.7)
#' 
#' 
#' # Combine several diagrams using grid graphics
#' vplayout = function(x, y){
#' 	return(viewport(layout.pos.row = x, layout.pos.col = y))
#' }
#' grid.newpage()
#' pushViewport(viewport(layout = grid.layout(ncol = 2, nrow = 2)))
#' for(i in 1:2){
#' 	for(j in 1:2){
#' 		pushViewport(vplayout(i, j))
#' 		bvenn(list(Set1 = sample(letters, 16), Set2 = sample(letters, 3+ 3*j), Set3 = sample(letters, 7)), add = TRUE, fontsize = 10)
#' 		upViewport()
#' 	}
#' }
#' 
#' 
#' @export
bvenn = function(sets, percentage = FALSE, colors = NULL, fontsize = 15, main = "", scale = 1, add = FALSE){
	n = length(sets)
	if(n < 2){
		stop("Please provide at least 2 lists")
	}
	
	if(!add){
		grid.newpage()
	}
	
	if(is.null(colors)){
		if(n == 2){
			colors = c("#377EB8", "#F4F422",  "#4DAF4A")
		}
		if(n > 2){
			colors = c("#E41A1C", "#377EB8", "#F4F422", "#984EA3", "#FF7F00", "#4DAF4A", "#A65628")
		}
	}
	
	if(n == 2 & length(colors) != 3){
		stop("Number of colors for 2 set diagram has to be 3")
	}
	if(n > 2 & length(colors) != 7){
		stop("Number of colors for 3 set diagram has to be 7")
	}
	
	if(main != ""){
		height = unit.c(unit(1.2, "lines"), unit(1, "npc") - unit(1.2, "lines"))
	}
	else{
		height = unit.c(unit(0, "lines"), unit(1, "npc"))
	}
	
	pushViewport(viewport(layout = grid.layout(ncol = 1, nrow = 2, heights = height), gp = gpar(fontsize = fontsize)))
	
	if(main != ""){
		pushViewport(vplayout(1, 1))
		grid.text(main, y = 0, vjust = 0, gp = gpar(fontface = "bold", cex = 1.3))
		popViewport()
	}
	
	pushViewport(vplayout(2, 1))
	if(n == 2){
		bubbleplot2l(sets, percentage = percentage, colors = colors, scale = scale)
	}
	if(n > 2){
		bubbleplot3l(sets[1:3], percentage = percentage, colors = colors, scale = scale)
	}
	popViewport()
	popViewport()
}







