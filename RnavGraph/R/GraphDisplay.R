## write NG_graph objcet onto the canvas
.graph2canvas <- function(ngEnv) {
	
	tkdelete(ngEnv$canvas, 'all')

		## resize graph
	width <- .tcl2num(tcl('winfo','width',ngEnv$canvas))
	height <- .tcl2num(tcl('winfo','height',ngEnv$canvas))
	
	ldist <- ngEnv$settings@interaction@labelDistRadius
	twoB <- 2*ngEnv$graph@border
	ngEnv$graph@xNodes <- (ngEnv$graph@xNodes - ngEnv$graph@border)/
			(ngEnv$graph@bbox[1]-twoB)*(width-twoB)+ngEnv$graph@border
	ngEnv$graph@yNodes <- (ngEnv$graph@yNodes-ngEnv$graph@border)/
			(ngEnv$graph@bbox[2]-twoB)*(height-twoB)+ngEnv$graph@border
	
	ngEnv$graph@xLabels = (ngEnv$graph@xLabels - ngEnv$graph@border) /
			(ngEnv$graph@bbox[1] - twoB) * (width - twoB)+ngEnv$graph@border
	ngEnv$graph@yLabels = (ngEnv$graph@yLabels - ngEnv$graph@border) /
			(ngEnv$graph@bbox[2]-twoB) * (height - twoB) + ngEnv$graph@border
	
	ngEnv$graph@bbox <- c(width,height)
	
	
	## update labels
	dx <- ngEnv$graph@xLabels - ngEnv$graph@xNodes 
	dy <- ngEnv$graph@yLabels - ngEnv$graph@yNodes
	
	dist <- apply(cbind(dx,dy),1,function(row)sqrt(sum(row^2)))
	
	ratio <- ldist/dist
	ngEnv$graph@xLabels <- ngEnv$graph@xNodes + ratio*dx
	ngEnv$graph@yLabels <- ngEnv$graph@yNodes + ratio*dy
	
	assign(ngEnv$cGraph,ngEnv$graph, envir = ngEnv)
	
#	.normalState(ngEnv)
#	cat("HERE\n")
	
	## edges
	
	## this can be programmed more efficient with a from-to-matrix
	## edges(ngEnv$graph@graph)
	
	active <- ngEnv$settings@color@edgeActive 
	width <- ngEnv$settings@display@lineWidth
	
	
	edgeM <- edgeMatrix(as(ngEnv$graph@graph,"graphNEL"))
	node <- nodes(ngEnv$graph@graph)
	
	if(dim(edgeM)[2]>0){
		for(i in 1:dim(edgeM)[2]) {
			if (edgeData(ngEnv$graph@graph,"visited",from=node[edgeM[1,i]], to=node[edgeM[2,i]])==TRUE) {
				tag <- c("graph", "edge", node[edgeM[1,i]], node[edgeM[2,i]],'visited')
				col <- ngEnv$settings@color@visitedEdge
			} else {	
				tag <- c("graph", "edge", node[edgeM[1,i]], node[edgeM[2,i]])
				col <- ngEnv$settings@color@notVisitedEdge
			}
			tkcreate(ngEnv$canvas, 'line',
					ngEnv$graph@xNodes[edgeM[1,i]],ngEnv$graph@yNodes[edgeM[1,i]],
					ngEnv$graph@xNodes[edgeM[2,i]],ngEnv$graph@yNodes[edgeM[2,i]],
					fill = col,
					activefill = active,
					width = width,
					tag = tag)
		}
	}	
	
	
	## draw nodes
	r <- ngEnv$settings@display@nodeRadius
	col <- ngEnv$settings@color@nodes
	active <- ngEnv$settings@color@nodesActive
	mapply(nodes(ngEnv$graph@graph),ngEnv$graph@xNodes,ngEnv$graph@yNodes,
			FUN = function(node,x,y) {
				tkcreate(ngEnv$canvas, 'oval', x-r, y-r, x+r, y+r,
						fill = col,
						activefill = active,
						tag = c("graph", "node", node)
				)
			})
	
	
	## draw labels
	col <-  ngEnv$settings@color@labels
	active <- ngEnv$settings@color@labelsActive 
	mapply(nodes(ngEnv$graph@graph),ngEnv$graph@xLabels,ngEnv$graph@yLabels,
			FUN = function(node,x,y) {
				tkcreate(ngEnv$canvas, 'text', x, y,
						text = node,
						fill = col,
						activefill = active,
						tag = c("graph", "label", node)
				)
			})
	
	
	
	## draw bullet
	r <- ngEnv$settings@display@bulletRadius
	col <- ngEnv$settings@color@bullet
	active <- ngEnv$settings@color@bulletActive
	coord <- .tcl2xy(tkcoords(ngEnv$canvas, paste('node &&', ngEnv$bulletState$from)))
	
	if(ngEnv$bulletState$percentage != 0) {
		to <- .tcl2xy(tkcoords(ngEnv$canvas, paste('node &&', ngEnv$bulletState$to)))
		coord <- coord + (coord-to)*ngEnv$bulletState$percentage
	}
	
	tkcreate(ngEnv$canvas, 'oval', coord[1]-r, coord[2]-r, coord[1]+r, coord[2]+r,
			fill = col,
			activefill = active,
			tag = "bullet")
	
	## create percentage text
	tkcreate(ngEnv$canvas, 'text', 20, 20, 
			text = floor(ngEnv$bulletState$percentage*100),
			fill = ngEnv$settings@color@labels,
			tag = "percentage")
	
}



## update NG_graph object from canvas state
.canvas2graph <- function(ngEnv) {
	## create a new graph object
	nodeMat <- matrix(.tcl2str(tcl('getNodeNames',ngEnv$canvas)), ncol = 5, byrow = FALSE)

	
	nodes <- nodeMat[,1]
	ngEnv$graph@xNodes <- as.numeric(nodeMat[,2])
	ngEnv$graph@yNodes <- as.numeric(nodeMat[,3])
	ngEnv$graph@xLabels <- as.numeric(nodeMat[,4])
	ngEnv$graph@yLabels <- as.numeric(nodeMat[,5])
	
	w <- .tcl2num(tcl('winfo','width',ngEnv$canvas))
	h <- .tcl2num(tcl('winfo','height',ngEnv$canvas))
	ngEnv$graph@bbox <- c(w,h)
	
	
	## edges
	ftL <- .tcl2str(tcl('getEdgeFromToList',ngEnv$canvas))
	
	
	if(length(ftL) == 0) {
		## no edges
		ngEnv$graph@graph <- as(newgraph(nodes, mat = matrix(rep(0,length(nodes)^2), ncol = length(nodes)), isAdjacency = TRUE ),"graphNEL")
		edgeDataDefaults(ngEnv$graph@graph,"visited") <- FALSE
	}else {
		## edges
		matFT <- matrix(ftL, ncol = 3, byrow = FALSE)
	
		ngEnv$graph@graph <- newgraph(nodes, mat = matFT[,1:2], isAdjacency = FALSE)
		edgeDataDefaults(ngEnv$graph@graph,"visited") <- FALSE
		apply(matFT, 1, FUN = function(row){
					edgeData(ngEnv$graph@graph,"visited", from = row[1], to = row[2]) <- (row[3] != '-1')
				})
	}
	
}



.edgeSeen <- function(ngEnv) {
	tag <- paste('edge && ', ngEnv$bulletState$from, ' && ', ngEnv$bulletState$to)
	tkitemconfigure(ngEnv$canvas, tag, fill = ngEnv$settings@color@visitedEdge)
	tcl(ngEnv$canvas,'addtag', 'visited', 'withtag', tag)
}


.arriveAtNode <- function(ngEnv) {
	ngEnv$bulletState$to <- ""
	.highlightAdj(ngEnv)
}


.leaveNode <- function(ngEnv) {
	.highlightCurrent(ngEnv)
}

.highlightAdj <- function(ngEnv, from = NULL) {
	if(is.null(from)) {
		from <- ngEnv$bulletState$from
	}
	
	adjNodeV <- adj(ngEnv$graph@graph, from)[[from]]
	if(length(adjNodeV)>0) {
		## highlight adjacent edges
		tkitemconfigure(ngEnv$canvas,'edge', width = ngEnv$settings@display@lineWidth)
		tkitemconfigure(ngEnv$canvas,paste('edge && ', from), width = ngEnv$settings@display@highlightedLineWidth)
		tkitemraise(ngEnv$canvas,paste('edge && ', from),'edge')
		
		## color adjacent nodes
		tkitemconfigure(ngEnv$canvas,'node', fill = ngEnv$settings@color@nodes) #, activefill = ngEnv$settings@color@nodesActive)
		
		adjNodes <- paste(adjNodeV, collapse = ' || ')
		
		tkitemconfigure(ngEnv$canvas,paste('node && (', adjNodes ,')'), fill = ngEnv$settings@color@adjNodes) #, activefill = ngEnv$settings@color@adjNodesActive)
		
		
		## labels
		tkitemconfigure(ngEnv$canvas,'label', fill = ngEnv$settings@color@labels) #, activefill = ngEnv$settings@color@labelsActive)
		tkitemconfigure(ngEnv$canvas,paste('label && (', adjNodes ,')'), fill = ngEnv$settings@color@adjLabels) #, activefill = ngEnv$settings@color@adjLabelsActive)
	}
}

.highlightCurrent <- function(ngEnv) {
	## edges
	tkitemconfigure(ngEnv$canvas,'edge', width = ngEnv$settings@display@lineWidth)
	tkitemconfigure(ngEnv$canvas,paste('edge && ', ngEnv$bulletState$from, ' && ', ngEnv$bulletState$to), width = ngEnv$settings@display@highlightedLineWidth)
	tkitemraise(ngEnv$canvas,paste('edge && ', ngEnv$bulletState$from, ' && ', ngEnv$bulletState$to), 'edge')
	
	## color adjacent nodes
	tkitemconfigure(ngEnv$canvas,'node', fill = ngEnv$settings@color@nodes) #, activefill = ngEnv$settings@color@nodesActive)
	tkitemconfigure(ngEnv$canvas, paste('node && (',ngEnv$bulletState$from, ' || ', ngEnv$bulletState$to,')'), fill =  ngEnv$settings@color@adjNodes) #, activefill = ngEnv$settings@color@adjLabelsActive)
	
	## labels
	tkitemconfigure(ngEnv$canvas,'label', fill = ngEnv$settings@color@nodes, activefill = ngEnv$settings@color@nodesActive)
	tkitemconfigure(ngEnv$canvas, paste('label && (',ngEnv$bulletState$from, ' || ', ngEnv$bulletState$to,')'), fill = ngEnv$settings@color@adjLabels) #, activefill = ngEnv$settings@color@adjLabelsActive)
	
	## so that any-leave does not recolor it back...
	ngEnv$tempColor <- ngEnv$settings@color@adjLabels
}

.normalState <- function(ngEnv) {
	tkitemconfigure(ngEnv$canvas, 'bullet', fill = ngEnv$settings@color@bullet) #, activefill = ngEnv$settings@color@bulletActive)
	
	tkitemconfigure(ngEnv$canvas, 'edge', fill = ngEnv$settings@color@notVisitedEdge)
	tkitemconfigure(ngEnv$canvas, 'visited', fill = ngEnv$settings@color@visitedEdge)
	
	if((ngEnv$bulletState$percentage == 0) && (ngEnv$bulletState$to == '')) {
		.highlightAdj(ngEnv)
	}else {
		.highlightCurrent(ngEnv)
	}
}


.highlightPath <- function(ngEnv,from,to) {
	tkitemconfigure(ngEnv$canvas,'edge', width = ngEnv$settings@display@lineWidth)
	tkitemconfigure(ngEnv$canvas, paste('edge && ', from, ' && ',to ,sep = ''), fill = ngEnv$settings@color@path, width=ngEnv$settings@display@lineWidth)
	.highlightAdj(ngEnv,to)
	tkitemconfigure(ngEnv$canvas, 'bullet', fill = ngEnv$settings@color@bullet)
	
	if(ngEnv$bulletState$from %in% c(from,to)) {
		tkitemconfigure(ngEnv$canvas, 'bullet', fill = ngEnv$settings@color@path)			
	}
	
	tkitemconfigure(ngEnv$canvas, paste('(label && ', to, ') || (node && ',to, ')'), fill = ngEnv$settings@color@path)
	tkitemconfigure(ngEnv$canvas, paste('node && ',from), fill = ngEnv$settings@color@path)
	
	## so that any-leave does not recolor it back...
	ngEnv$tempColor <- ngEnv$settings@color@path
}



.walkEdge <- function(ngEnv, backward = FALSE, path = FALSE) {
	## This function should eventually use tcl multi threading
	## Note: it looks at the percentage
	
	ds <- 1/ngEnv$settings@interaction@NSteps ## step size
	
	if(backward) {
		nStepsToGo <- floor(ngEnv$settings@interaction@NSteps*ngEnv$bulletState$percentage)
		pSteps <- seq(ngEnv$bulletState$percentage-ds, 0, length.out = nStepsToGo)
	}else {
		nStepsToGo <- floor(ngEnv$settings@interaction@NSteps*(1-ngEnv$bulletState$percentage))
		pSteps <- seq(ngEnv$bulletState$percentage+ds, 1, length.out = nStepsToGo)
	}
	
	dt <- ngEnv$settings@interaction@animationTime/ngEnv$settings@interaction@NSteps
	
	
	xyfrom <- .tcl2xy(tkcoords(ngEnv$canvas,
					paste(ngEnv$bulletState$from,'&& node')))
	xyto <- .tcl2xy(tkcoords(ngEnv$canvas,
					paste(ngEnv$bulletState$to,'&& node')))
	dNodeXY <- xyto-xyfrom
	
	for(p in pSteps) {
		xy <- xyfrom + dNodeXY*p
		xybullet <- .tcl2xy(tkcoords(ngEnv$canvas,'bullet'))
		
		dMoveX <- xy[1]-xybullet[1]
		dMoveY <- xy[2]-xybullet[2]
		tkmove(ngEnv$canvas,'bullet', dMoveX, dMoveY)
		
		
		if(!(p %in% c(0,1))){
			ngEnv$bulletState$percentage <- p
			.updatePlots(ngEnv)
			tcl("update", "idletasks")
			Sys.sleep(dt)
		}
		tkitemconfigure(ngEnv$canvas, 'percentage', text = floor((p*100)%%100))
	}
	
	## update state
	ngEnv$bulletState$percentage <- 0
	.edgeSeen(ngEnv)
	
	if(!backward){
		ngEnv$bulletState$from <- ngEnv$bulletState$to
	}
	ngEnv$bulletState$to = ""	
	if(!path){
		.arriveAtNode(ngEnv)	
	}
	.updatePlots(ngEnv)
	tcl("update", "idletasks")
}


.moveBullet <- function(ngEnv,x,y) {
	x <- as.numeric(x)
	y <- as.numeric(y)
	
	## mouse change
	dx <- x - ngEnv$mouseHistory$x
	dy <- y - ngEnv$mouseHistory$y
	
	## effective moving difference of bullet
	if((ngEnv$bulletState$percentage == 0) && (ngEnv$bulletState$to == '')) {
		if(ngEnv$mouseHistory$Click) {
			adjNodes <- adj(ngEnv$graph@graph,ngEnv$bulletState$from)[[ngEnv$bulletState$from]]
			
			if(length(adjNodes) > 0) {
				## select node
				dxClick <- x - ngEnv$mouseHistory$xClick
				dyClick <- y - ngEnv$mouseHistory$yClick
				
				if(sqrt(dxClick^2+dyClick^2)>ngEnv$settings@interaction@dragSelectRadius) {
					## choose direction
					ii <- which.min(abs(.tcl2num(tcl('minThetaNode',ngEnv$canvas, ngEnv$bulletState$from, paste(adjNodes,collapse = ' '), dxClick, dyClick))))
					
					ngEnv$bulletState$to <- adjNodes[ii]
					.leaveNode(ngEnv)
					
					## move bullet to the right place
					xyfrom <- .tcl2xy(tkcoords(ngEnv$canvas,
									paste(ngEnv$bulletState$from,'&& node')))
					xyto <- .tcl2xy(tkcoords(ngEnv$canvas,
									paste(ngEnv$bulletState$to,'&& node')))
					dNodeXY <- xyto-xyfrom
					
					## find percentage
					if(which.max(abs(dNodeXY)) == 1) {
						## when dNodeX is bigger use x-direction sensitivity
						ngEnv$bulletState$percentage <- (x-xyfrom[1])/dNodeXY[1]
					}else {
						## when dy is bigger use y-direction sensitivity
						ngEnv$bulletState$percentage <- (y-xyfrom[2])/dNodeXY[2]
					}
					
					if(ngEnv$bulletState$percentage >= 1) {
						## arrive at new node
						
						ngEnv$bulletState$percentage <- 0
						.edgeSeen(ngEnv)
						ngEnv$bulletState$from <- ngEnv$bulletState$to
						.arriveAtNode(ngEnv)
						
						xyfrom <- xyto
						xy <- xyfrom + (xyto-xyfrom)*ngEnv$bulletState$percentage
						xybullet <- .tcl2xy(tkcoords(ngEnv$canvas,'bullet'))
						
						dMoveX <- xy[1]-xybullet[1]
						dMoveY <- xy[2]-xybullet[2]
						ngEnv$mouseHistory$Click <- FALSE
						
					}else if(ngEnv$bulletState$percentage <= 0) {
						ngEnv$bulletState$percentage <- 0
						.arriveAtNode(ngEnv)
						
						xy <- xyfrom + (xyto-xyfrom)*ngEnv$bulletState$percentage
						xybullet <- .tcl2xy(tkcoords(ngEnv$canvas,'bullet'))
						
						dMoveX <- xy[1]-xybullet[1]
						dMoveY <- xy[2]-xybullet[2]
						ngEnv$mouseHistory$Click <- FALSE
					}else {
						
						xybullet <- .tcl2xy(tkcoords(ngEnv$canvas,'bullet'))
						
						dxyMove <- xyfrom + dNodeXY* ngEnv$bulletState$percentage - xybullet
						
						dMoveX <- dxyMove[1]
						dMoveY <- dxyMove[2]
					}
					
				}else {
					dMoveX <- dx
					dMoveY <- dy 
				} ## sqrt(dist>sel...)
			}else {
				## no adjacent
				dMoveX <- 0
				dMoveY <- 0      
			}
		}else { ## no wait until click
			dMoveX <- 0
			dMoveY <- 0     
		}
	}else { ## percentage != 0
		## traversal on edge
		xyfrom <- .tcl2xy(tkcoords(ngEnv$canvas,
						paste(ngEnv$bulletState$from,'&& node')))
		xyto <- .tcl2xy(tkcoords(ngEnv$canvas,
						paste(ngEnv$bulletState$to,'&& node')))
		dNodeXY <- xyto-xyfrom
		
		## find percentage change
		if(which.max(abs(dNodeXY)) == 1) {
			## when dNodeX is bigger use x-direction sensitivity
			dpercentage <- dx/dNodeXY[1]
		}else {
			## when dy is bigger use y-direction sensitivity
			dpercentage <- dy/dNodeXY[2]
		}
		
		ngEnv$bulletState$percentage <- ngEnv$bulletState$percentage + dpercentage
		p <- ngEnv$bulletState$percentage
		if(p >0 && p<1) {
			## traversal on edge
			
			## change nothing
			
		}else if(p<=0) {
			## go back to from node
			ngEnv$bulletState$percentage <- 0
			.arriveAtNode(ngEnv)
			
		}else if(p>=1) {
			## arrive at new node
			ngEnv$bulletState$percentage <- 0
			.edgeSeen(ngEnv)
			ngEnv$bulletState$from <- ngEnv$bulletState$to
			.arriveAtNode(ngEnv)
			xyfrom <- xyto
			
		}else {
			## not possible: reset
			ngEnv$bulletState$percentage <- 0
		}
		
		xy <- xyfrom + (xyto-xyfrom)*ngEnv$bulletState$percentage
		xybullet <- .tcl2xy(tkcoords(ngEnv$canvas,'bullet'))
		
		dMoveX <- xy[1]-xybullet[1]
		dMoveY <- xy[2]-xybullet[2]
		
	}
	
	if(any(dMoveX != 0 , dMoveY != 0)) {
		tkmove(ngEnv$canvas, 'bullet', dMoveX, dMoveY)
		tcl('update', 'idletasks')
		.updatePlots(ngEnv)
	}
	tkitemconfigure(ngEnv$canvas,'percentage', text = floor(ngEnv$bulletState$percentage*100))
	
	ngEnv$mouseHistory$x <- x
	ngEnv$mouseHistory$y <- y
	
	if(ngEnv$bulletState$percentage != 0){
		ngEnv$mouseHistory$Click <- FALSE
	}
}

.moveNode <- function(ngEnv,x,y) {
	x <- as.numeric(x)
	y <- as.numeric(y)
	dx <- x - ngEnv$mouseHistory$x
	dy <- y - ngEnv$mouseHistory$y
	
	## get node name
	node <- .tcl2str(tkgettags(ngEnv$canvas, 'current'))[3]
	
	## move Node
	tkmove(ngEnv$canvas, 'current', dx, dy)
	## move Label
	tkmove(ngEnv$canvas, paste("label && ", node ), dx, dy)
	
	
	## move Edges
	adjNodesV <- adj(ngEnv$graph@graph,node)[[node]]
	adjNodes <- paste(adjNodesV, collapse = ' ')
	
	tcl('adjustEdges',ngEnv$canvas,node,adjNodes)
	
	if((ngEnv$bulletState$from == node && ngEnv$bulletState$to %in% adjNodesV) ||
			(ngEnv$bulletState$to == node && ngEnv$bulletState$from %in% adjNodesV)) {
		xy1 <- .tcl2xy(tkcoords(ngEnv$canvas,paste('node &&', ngEnv$bulletState$from)))
		xy2 <- .tcl2xy(tkcoords(ngEnv$canvas,paste('node &&', ngEnv$bulletState$to)))
		
		xy <- xy1 + (xy2-xy1)*ngEnv$bulletState$percentage
		bxy <- .tcl2xy(tkcoords(ngEnv$canvas,'bullet'))
		dxy <- xy-bxy
		
		if(any(dxy[1] != 0 , dxy[2] != 0)) {
			tkmove(ngEnv$canvas,'bullet',dxy[1],dxy[2])
			.updatePlots(ngEnv)
		}	
	}
	
	## update mouse history
	ngEnv$mouseHistory$x <- x
	ngEnv$mouseHistory$y <- y
	
	ngEnv$mouseHistory$Click <- FALSE
}

.moveLabel <- function(ngEnv,x,y) {
	x <- as.numeric(x)
	y <- as.numeric(y)
	
	## get label name
	label <- .tcl2str(tkgettags(ngEnv$canvas, 'current'))[3]
	nodexy <- .tcl2xy(tkcoords(ngEnv$canvas, paste('node && ',label)))
	
	vec <- c(x,y)-nodexy
	len <- sqrt(sum(vec^2))
	if(len <= ngEnv$settings@interaction@labelDistRadius) {
		tkcoords(ngEnv$canvas,'current',x,y)
	} else {
		tkcoords(ngEnv$canvas,'current',
				nodexy[1]+vec[1]/len*ngEnv$settings@interaction@labelDistRadius,
				nodexy[2]+vec[2]/len*ngEnv$settings@interaction@labelDistRadius)
	}
}