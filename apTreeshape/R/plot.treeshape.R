"plot.treeshape" <-
function(x, y, ...){

plotone <- function(tree, ...){
	n <- nrow(tree$merge)+400
	hc <-  hclust(d = dist(runif(n), method = "euclidean"), method = "ward")
	hc$merge <- tree$merge
	
	hc$height <- 1:nrow(tree$merge)
	
	descendant <- smaller.clade.spectrum(tree)[,1]
	current <- 1
	for (i in 2:(nrow(tree$merge)+1)) {
		if (sum(descendant==i)!=0) {
			descendant[descendant==i] <- current
			current <- current+1
		}
	}
	for (i in 1:length(descendant)) {
		hc$height[length(descendant)-i+1]=descendant[i]
	}

	#return(hc$height)
	
	hc$labels <- tree$names
	hc <- as.dendrogram(hc)
	
	mar<-par()$mar
	par(mar=c(1,1,1,10))
	plot(hc, horiz=TRUE, axes=FALSE, ...)
	
	par(mar=mar)
}
	tree1<-x
	
	if (missing(y)){
		plotone(tree1)
		#text(x=-2, y=1, label=tree1$names[1])
	}
	else{
		tree2<-y
		if (class(tree2)!='treeshape') {
			stop("invalid arguments")
		}
	
		layout(t(c(1,2)))
		plotone(tree1, ...)
		plotone(tree2, ...)
	layout(1)
	}
}

