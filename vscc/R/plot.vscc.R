plot.vscc <-
function(x, ...){
	classcolours <- rainbow(length(unique(x$bestmod$class)))
	pairs(x$top, col=classcolours[x$best$class])
}
