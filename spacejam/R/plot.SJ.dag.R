plot.SJ.dag <-
function(x, which = NULL, layout=NULL,...){
	if(is.null(which)){
		which <- which.min(x$bic)
	}
	if(is.null(layout)){
		set.seed(20)
		layout <- layout.fruchterman.reingold(x$graph[[which]])
	}
	plot(x$graph[[which]],layout=layout, edge.color = "gray50", vertex.color = "red", vertex.size = 3, vertex.label = NA, edge.arrow.size=0.4, ...)
	invisible(layout)
}
