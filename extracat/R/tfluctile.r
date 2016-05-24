
drawt = function(tree, vertical = FALSE, lwd = 1.2,line.col= 1, vp = NULL, ... ){
	if(is.null(tree)){
		return(invisible(TRUE))	
	}
	ng <- length(tree$order)
	maxht <- max(tree$height)
	nht <- nodeheight( tree$merge, tree$order, tree$height/maxht	)
	#nrg <- noderange( tree$merge, tree$order )
	#nleft <- lapply(nrg, function(z) mean(z[[1]])/ng )
	#nright <- lapply(nrg, function(z) mean(z[[2]])/ng )
	nxc <- nodexc(tree$merge, tree$order)
	if( !vertical ){
	mapply( function(x,hlr, h){
		grid.lines( x = c(x[1],x[2]), y = c(h,h), vp = vp, gp = gpar(lwd = lwd, col = line.col) )
		grid.lines( x = c(x[1],x[1]), y = c(h-hlr[[1]],h), vp = vp , gp = gpar(lwd = lwd, col = line.col))
		grid.lines( x = c(x[2],x[2]), y = c(h-hlr[[2]],h) , vp = vp, gp = gpar(lwd = lwd, col = line.col))
		}, x = nxc , hlr = nht, h = as.list(tree$height/maxht))
	}else{
		mapply( function(x,hlr, h){
		grid.lines( y = 1-c(x[1],x[2]), x = 1-c(h,h) , vp = vp, gp = gpar(lwd = lwd, col = line.col))
		grid.lines( y = 1-c(x[1],x[1]), x = 1-c(h-hlr[[1]],h) , vp = vp, gp = gpar(lwd = lwd, col = line.col))
		grid.lines( y = 1-c(x[2],x[2]), x = 1-c(h-hlr[[2]],h) , vp = vp, gp = gpar(lwd = lwd, col = line.col))
		}, x = nxc , hlr = nht, h = as.list(tree$height/maxht))
	}
invisible(TRUE)
	
}


tfluctile = function(x, tree = NULL, dims = c(1,2), tw = 0.2, border = NULL, shape = "r", dir = "b", just = "c",  tile.col= hsv(0.1,0.1,0.1,alpha=0.6),
 bg.col = "lightgrey", vp = NULL, lab.opt = list(), ... ){
	stopifnot( !is.null( attr(x,"tree") ) | !is.null(tree) )
	
	
	
	
	if("lwd" %in% names(lab.opt)){
		lwd <- lab.opt$lwd
	}else{
		lwd <- 1.2		
	}
	if("line.col" %in% names(lab.opt)){
		line.col <- lab.opt$line.col
	}else{
		line.col <- 1		
	}
	
	if(is.null(border)){
		border <- 0.05	
	}
	if( !is.null(tree) ){
		tree1 <- tree[[ dims[1] ]]	
		tree2 <- tree[[ dims[2] ]]	
	}else{
		tree1 <- attr(x,"tree")[[ dims[1] ]]	
		tree2 <- attr(x,"tree")[[ dims[2] ]]
	}
	if(length(dim(x)) > 2){
		x <- as.data.frame(as.table(x))
		x <- xtabs(x$Freq~x[,dims[1]]+x[,dims[2]])
	}
	
	if( !all(is.na(tree1)) ){
		
		m1 <- match(tree1$label,rownames(x))
		tree1$order <- m1[tree1$order]
		w1 <- which(tree1$merge < 0)
		tree1$merge[w1] <- -m1[-tree1$merge[w1]]
		tree1$label <- rownames(x)
		
		x <- x[tree1$order,]
	}
	if( !all(is.na(tree2)) ){
		m2 <- match(tree2$label,colnames(x))
		tree2$order <- m2[tree2$order]
		w2 <- which(tree2$merge < 0)
		tree2$merge[w2] <- -m2[-tree2$merge[w2]]
		tree2$label <- colnames(x)
		
		x <- x[,tree2$order]
	}
	
	
	lefts <- ifelse( all(is.na(tree1)), 0, tw )
	tops <- ifelse( all(is.na(tree2)), 0, tw )
	
	
	
	
	
	if(is.null(vp)){
		grid.newpage()	
	}else{
		pushViewport(vp)	
	}
	border <- border/2
	vp00 <- viewport(x=lefts,y=1-tops, just=c("left","top"), width = 1-lefts, height = 1-tops)
	
	if(!all(is.na(tree1))){
	vp0L <- viewport(x=0,y=1-tops, just=c("left","top"), width = lefts, height = (1-tops))
	pushViewport(vp0L)
	pushViewport( viewport(0.5,0.5,0.96,1-2*border) )
	drawt(tree1,vertical=TRUE, lwd = lwd, line.col = line.col, vp = NULL)
	upViewport(2)
	}
	if(!all(is.na(tree2))){
	vp0T <- viewport(x=lefts,y=1, just=c("left","top"), width = (1-lefts), height = tops)	
	pushViewport(vp0T)
	pushViewport( viewport(0.5,0.5,1-2*border,0.96) )
	drawt(tree2,vertical=FALSE, lwd = lwd, line.col = line.col, vp = NULL)
	upViewport(2)
	}
	pushViewport(vp00)
	fluctile(x, add = TRUE, tile.col = tile.col, bg.col = bg.col, lab.opt = lab.opt, shape = shape, dir = dir, just = just, border = 2*border)
	upViewport()
		
	return(invisible(TRUE))
}

