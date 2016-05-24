graphVAR1 <- function(sparseA, sparseP, type="TSCG", side="left", prune=TRUE, nNames=NULL, main=NULL, vertex.color.T0="lightcyan2", vertex.color.T1="lightcyan2", vertex.frame.color="steelblue", vertex.label.cex=-1, vertex.label.color.T0="black", vertex.label.color.T1="black", vertex.label.font=1.5, vertex.size=-1, edge.arrow.size=-1, edge.width=-1, ...){

	############################################################################################################
	#
	# DESCRIPTION:
	# Plot temporal relations as implied by A (the matrix with regression coefficient of the VAR(1)) model.
	#
	# ARGUMENTS:
	#
	# -> sparseA                 : Matrix A of regression parameters, which is assumed to be sparse.
	# -> sparseP                 : Matrix P of precision of the error, which is assumed to be sparse.
	# -> type                    : A 'character' indicating what should be plotted. If type='TSCG'. the time
	#                              series chain graph is plotted, while type type='Aonly' limits this graph
	#                              to the temporal relations. If type='globalPC' or type='contempPC', the 
	#                              global or contemporaneous (respectively) partial correlation graph is
	#                              plotted.
	# -> side                    : A 'character' indicating wether the contemporaneous dependencies should be plotted
	#                              on the 'left' (time t) or the 'right' (time t+1) side. Only active when type='TSCG'.
	# -> prune                   : A 'logical' indicating whether to remove covariates without any 
	#                              temporal relations (as implied by 'sparseA').
	# -> nNames                  : A 'character' containing the covariate names to written inside the nodes.
	# -> main                    : The 'character' to be plotted as title above the graph.
	# -> vertex.color.T0         : Color of nodes at time point t.
	# -> vertex.color.T1         : Color of nodes at time point t+1. Ignored when type='globalPC' or type='contempPC'.
	# -> vertex.frame.color      : Refer to 'plot.igraph'.
	# -> vertex.label.cex        : Refer to 'plot.igraph'.
	# -> vertex.label.color.T0   : Color of the node label at time point t.
	# -> vertex.label.color.T1   : Color of the node label at time point t+1. Ignored when type='globalPC' or type='contempPC'.
	# -> vertex.label.font       : Refer to 'plot.igraph'.
	# -> vertex.size             : Refer to 'plot.igraph'.
	# -> edge.arrow.size         : Refer to 'plot.igraph'.
	# -> edge.width              : Refer to 'plot.igraph'.
	# -> ...                     : Other arguments to be passed to 'plot.igraph'.
	# 
	# DEPENDENCIES:
	# library(igraph)	    # functions: delete.edges, E, ecount, plot.igraph, graph.adjacency, maximum.cardinality.search
	#
	# NOTES:
   	# - igraph does not support visualization of mixed graphs, including both directed and undirected edges.
	# as a consequence the time-series chain graph is undirected: the edges from t to t+1 should be thought of as directed. 
	# - if "white" labels, color combination "navy" and "blue" looks nice. 
	# - if "black" labels, color combination "tomato1" and "red" looks nice. 
	#     	
	############################################################################################################

	# input check
	if (as.character(class(sparseA)) != "matrix"){ stop("Input (sparseA) is of wrong class.") }
	if (nrow(sparseA) != ncol(sparseA)){ stop("Matrix sparseA is not square.") }
	if (as.character(class(sparseP)) != "matrix"){ stop("Input (sparseP) is of wrong class.") }
	if (nrow(sparseA) != ncol(sparseP)){ stop("Matrix sparseP is not square.") }
	if (as.character(class(prune)) != "logical"){ stop("Input (prune) is of wrong class.") }
	if (as.character(class(vertex.size)) != "numeric"){ stop("Input (vertex.size) is of wrong class.") }
	if (as.character(class(vertex.color.T0)) != "character"){ stop("Input (vertex.color.T0) is of wrong class.") }
	if (as.character(class(vertex.color.T1)) != "character"){ stop("Input (vertex.color.T1) is of wrong class.") }
	if (as.character(class(vertex.label.color.T0)) != "character"){ stop("Input (vertex.label.color.T0) is of wrong class.") }
	if (as.character(class(vertex.label.color.T1)) != "character"){ stop("Input (vertex.label.color.T1) is of wrong class.") }

	# if no covariate names are specified the columns and row names of A are given names 1, 2, et cetera
	if (is.null(nNames)){ nNames <- as.character(1:nrow(sparseA)) }

	if (type=="TSCG"){     
		# prune covariate without connections (according to sparseA and sparseP)
		if (prune){ 
			idRemoveA <- intersect(which(apply(sparseA, 1, function(Z){ all(Z == 0) })), which(apply(sparseA, 2, function(Z){ all(Z == 0) })))
			diag(sparseP) <- 0 
			idRemoveP <- intersect(which(apply(sparseP, 1, function(Z){ all(Z == 0) })), which(apply(sparseP, 2, function(Z){ all(Z == 0) })))
			idRemove <- intersect(idRemoveA, idRemoveP)
			if (length(idRemove) > 0){ 
				sparseA <- sparseA[-idRemove, -idRemove]
				sparseP <- sparseP[-idRemove, -idRemove] 
				nNames <- nNames[-idRemove]
			}
		}
		nNodes <- nrow(sparseA)
		if (vertex.label.cex <= 0){ vertex.label.cex <- max(6*(nrow(sparseA))^(-0.8), 0.1) }
		if (vertex.size <= 0){ vertex.size <- max(75*(nrow(sparseA))^(-0.7), 1) }
		if (is.null(main)){ main <- "VAR(1) time-series chain graph" }

		# reshuffle nodes
		CIG <- CIGofVAR1(sparseA, sparseP, type="global")        
		shuffle <- maximum.cardinality.search(graph.adjacency(CIG, "undirected"))$alpham1
		sparseA <- sparseA[shuffle, shuffle]
		sparseP <- sparseP[shuffle, shuffle]
		nNames <- nNames[shuffle]

		# rectangular plot layout
		grid <- rbind(cbind(-1, 1:nrow(sparseA)), cbind(1, 1:nrow(sparseA)))

		# contempory independencies
		contCIG <- sparseP 
		contCIG[sparseP != 0] <- 1
		diag(contCIG) <- 0
		contCIG[lower.tri(contCIG)] <- 0
        
		if (side == "right"){
			# adjacency matrix 
			adjMat <- rbind(cbind(0*sparseA, t(sparseA)), cbind(0*sparseA, contCIG))
			adjMat[adjMat != 0] <- 1
            
			# size of the edges       
			edge.widthA <- abs(t(sparseA)[which(t(sparseA) != 0, arr.ind=TRUE)]) 
			edge.widthP <- abs(sparseP[which(contCIG != 0, arr.ind=TRUE)])
			edge.widthA  <- edge.widthA  / max(edge.widthA)
			edge.widthP  <- edge.widthP  / max(edge.widthP)
			if (edge.width <= 0){ 
				edge.width <- 40 * (nrow(sparseA)^(-0.8)) * c(edge.widthA, edge.widthP)  
			} else { 
				edge.width <- edge.width * c(edge.widthA, edge.widthP)  
			}
        
			# convert adjacency matrix to graph-object
			gObj <- graph.adjacency(adjMat, mode="undirected")

			# arrow heads of curved edges should be zero width
			curved <- rep(0, sum(adjMat))
			curved[-c(1:sum(adjMat[1:nNodes,]))] <- -1        
			negEdges <- which(sparseA[which(sparseA != 0)] < 0)
			negEdges <- c(negEdges, length(which(sparseA != 0)) + which(sparseP[which(contCIG != 0, arr.ind=TRUE)] < 0))
			posEdges <- which(sparseA[which(sparseA != 0)] > 0)
			posEdges <- c(posEdges, length(which(sparseA != 0)) + which(sparseP[which(contCIG != 0, arr.ind=TRUE)] > 0))
			igraph::E(gObj)[negEdges]$style <- 2
			igraph::E(gObj)[posEdges]$style <- 1
    
			# make plot
			plot(gObj, edge.lty=0, edge.arrow.size=0, layout=grid, vertex.shape="none", vertex.color="white", main=main, margin=c(0, -0.5, 0.1, 0.3), vertex.label.color="white")
			plot(gObj, add=TRUE, layout=grid, vertex.size=rep(vertex.size, 2*nNodes), vertex.label.font=vertex.label.font, edge.width=edge.width,
				vertex.color=c(rep(vertex.color.T0, nNodes), rep(vertex.color.T1, nNodes)), vertex.frame.color=vertex.frame.color, 
				vertex.label=c(nNames, nNames), vertex.label.cex=vertex.label.cex, vertex.label.color=c(rep(vertex.label.color.T0, nNodes), 
				rep(vertex.label.color.T1, nNodes)),  edge.color="black", edge.width=edge.width, vertex.label.family="sans", 
				edge.curved=curved, edge.lty=igraph::E(gObj)$style, margin=c(0, -0.5, 0.1, 0.3), ...)
			text(-1, -1.2, "t", cex=1.2, font=3)
			text(1, -1.2, expression(paste(italic(t), "+1 ", sep="")), cex=1.2)
		}
        
		# make plot
		if (side == "left"){
			# adjacency matrix 
			adjMat <- rbind(cbind(contCIG, t(sparseA)), cbind(0*sparseA, 0*sparseA))
			adjMat[adjMat != 0] <- 1
            
			# size of the edges       
			edge.widthA <- abs(t(sparseA)[which(t(sparseA) != 0, arr.ind=TRUE)]) 
			edge.widthP <- abs(sparseP[which(contCIG != 0, arr.ind=TRUE)])
			edge.widthA  <- edge.widthA  / max(edge.widthA)
			edge.widthP  <- edge.widthP  / max(edge.widthP)            
			edge.widthBoth <- rep(NA, sum(adjMat))
			for (j in 1:nrow(sparseP)){
				slh <- which(adjMat[j, c(1:nrow(sparseP))]==1)     
				if (length(slh)){ edge.widthBoth[sum(!is.na(edge.widthBoth)) + 1:length(slh)] <- edge.widthP[1:length(slh)]; edge.widthP <- edge.widthP[-c(1:length(slh))] }
				slh <- which(adjMat[j, -c(1:nrow(sparseP))]==1)
				if (length(slh)){ edge.widthBoth[sum(!is.na(edge.widthBoth)) + 1:length(slh)] <- edge.widthA[1:length(slh)]; edge.widthA <- edge.widthA[-c(1:length(slh))] }
			}            
			if (edge.width <= 0){ 
				edge.width <- 40 * (nrow(sparseA)^(-0.8)) * edge.widthBoth  
			} else { 
				edge.width <- edge.width * edge.widthBoth  
			}
        
			# convert adjacency matrix to graph-object
			gObj <- graph.adjacency(adjMat, mode="undirected")

			# arrow heads of curved edges should be zero width
			curved <- rep(0, sum(adjMat))
			slh <- which(adjMat[1, c(1:nrow(sparseP))]==1)
			if (length(slh) > 0){ idCurved <- 1:length(slh) } else { idCurved <- numeric() }
			for (j in 2:nrow(sparseP)){
				slh <- which(adjMat[j, c(1:nrow(sparseP))]==1)
				if (length(slh) > 0){
					idCurved <- c(idCurved, 1:length(slh) + sum(adjMat[1:(j-1),]))
				}
			}
			curved[idCurved] <- 1

			negEdges <- numeric()
			posEdges <- numeric()
			for (j in 1:nrow(sparseP)){
				slh <- which(adjMat[j, c(1:nrow(sparseP))]==1)
				slh2 <- length(negEdges) + length(posEdges)
				if (length(slh) > 0){ negEdges <- c(negEdges, slh2 + which(sparseP[j, slh] < 0)); 
                                		        posEdges <- c(posEdges, slh2 + which(sparseP[j, slh] > 0)) }
				slh <- which(adjMat[j, -c(1:nrow(sparseP))]==1)
				slh2 <- length(negEdges) + length(posEdges)                
				if (length(slh) > 0){ negEdges <- c(negEdges, slh2 + which(t(sparseA)[j, slh] < 0)); 
                                		        posEdges <- c(posEdges, slh2 + which(t(sparseA)[j, slh] > 0)) }
			}            
			igraph::E(gObj)[negEdges]$style <- 2
			igraph::E(gObj)[posEdges]$style <- 1
       
			plot(gObj, edge.lty=0, edge.arrow.size=0, layout=grid, vertex.shape="none", vertex.color="white", main=main, margin=c(0, 0.2, 0.1, -0.1), vertex.label.color="white")
			plot(gObj, add=TRUE, layout=grid, vertex.size=rep(vertex.size, 2*nNodes), vertex.label.font=vertex.label.font, edge.width=edge.width,
				vertex.color=c(rep(vertex.color.T0, nNodes), rep(vertex.color.T1, nNodes)), vertex.frame.color=vertex.frame.color, 
				vertex.label=c(nNames, nNames), vertex.label.cex=vertex.label.cex, 
				vertex.label.color=c(rep(vertex.label.color.T0, nNodes), rep(vertex.label.color.T1, nNodes)),  edge.color="black",
				edge.width=edge.width, vertex.label.family="sans", edge.curved=curved, edge.lty=igraph::E(gObj)$style, margin=c(0, 0.2, 0.1, -0.1), ...)    
			text(-1, -1.2, "t", cex=1.2, font=3)
			text(1, -1.2, expression(paste(italic(t), "+1 ", sep="")), cex=1.2)
		}
	}
    
	if (type=="Aonly"){
		# prune covariate without connections (according to sparseA)
		if (prune){ 
			idRemove <- intersect(which(apply(sparseA, 1, function(Z){ all(Z == 0) })), which(apply(sparseA, 2, function(Z){ all(Z == 0) })))
			if (length(idRemove) > 0){ 
				sparseA <- sparseA[-idRemove, -idRemove]
				sparseP <- sparseP[-idRemove, -idRemove]
				nNames <- nNames[-idRemove]
			}
		}   
		nNodes <- nrow(sparseA)
		if (vertex.label.cex <=0 ){ vertex.label.cex <- max(6*(nrow(sparseA))^(-0.8), 0.1) }
		if (vertex.size <= 0){ vertex.size <- max(75*(nrow(sparseA))^(-0.7), 1) }
		if (is.null(main)){ main <- "Cross-time relations of VAR(1) model" }

		# reshuffle nodes
		CIG <- CIGofVAR1(sparseA, diag(diag(sparseP)), type="global")        
		shuffle <- maximum.cardinality.search(graph.adjacency(CIG, "undirected"))$alpham1
		sparseA <- sparseA[shuffle, shuffle]
		nNames <- nNames[shuffle]
		grid <- rbind(cbind(-1, 1:nrow(sparseA)), cbind(1, 1:nrow(sparseA)))
    
		# adjacency matrix 
		adjMat <- rbind(cbind(0*sparseA, t(sparseA)), cbind(0*sparseA, 0*t(sparseA)))
		adjMat[adjMat != 0] <- 1

		# convert adjacency matrix to graph-object
		gObj <- graph.adjacency(adjMat)

		# select pos and neg edges
		negEdges <- which(sparseA[which(sparseA != 0)] < 0)
		posEdges <- which(sparseA[which(sparseA != 0)] > 0)
		igraph::E(gObj)[negEdges]$style <- 2
		igraph::E(gObj)[posEdges]$style <- 1
           
		# size of the edges       
		edge.widthA <- abs(t(sparseA)[which(t(sparseA) != 0, arr.ind=TRUE)]) 
		edge.widthA  <- edge.widthA  / max(edge.widthA)
		if (edge.width <= 0){ 
			edge.width <- 40 * (nrow(sparseA)^(-0.8)) * edge.widthA  
		} else { 
			edge.width <- edge.width * edge.widthA  
		}
		if (edge.arrow.size <= 0){ 
			edge.arrow.size <- 10 * (nrow(sparseA)^(-0.8)) * edge.widthA 
		} else { 
			edge.arrow.size <- edge.arrow.size * edge.widthA  
		}

		# make plot
		plot(gObj, edge.lty=0, edge.arrow.size=0, layout=grid, vertex.shape="none", vertex.label=c(nNames, nNames), vertex.label.family="sans", main=main, vertex.label.cex=vertex.label.cex, ...)
		for (e in seq_len(ecount(gObj))){
			graph2 <- delete.edges(gObj, igraph::E(gObj)[(1:ecount(gObj))[-e]])
			plot(graph2, edge.arrow.size=edge.arrow.size[e], layout=grid, vertex.size=rep(vertex.size, 2*nNodes), vertex.label.font=vertex.label.font, 
				vertex.color=c(rep(vertex.color.T0, nNodes), rep(vertex.color.T1, nNodes)), vertex.frame.color=vertex.frame.color,
				layout=grid, edge.color="black", vertex.label.cex=vertex.label.cex, edge.width=edge.width[e], vertex.label.family="sans", edge.lty=igraph::E(graph2)$style,     
				vertex.label=c(nNames, nNames), add=TRUE, ...)        
		}        
		text(-1, -1.2, "t", cex=1.2, font=3)
		text(1, -1.2, expression(paste(italic(t), "+1 ", sep="")), cex=1.2)    
	}
    
	if (type=="globalPC"){     
		# get global CIG
		CIG <- CIGofVAR1(sparseA, sparseP, type="global")  
		diag(CIG) <- 0         

		# prune covariate without connections (according to sparseA and sparseP)
		if (prune){ 
			diag(CIG) <- 0         
			idRemove <- intersect(which(apply(CIG, 1, function(Z){ all(Z == 0) })), which(apply(CIG, 2, function(Z){ all(Z == 0) })))
			if (length(idRemove) > 0){ 
				sparseA <- sparseA[-idRemove, -idRemove]
				sparseP <- sparseP[-idRemove, -idRemove] 
				nNames <- nNames[-idRemove]
			}
		}
		nNodes <- nrow(sparseA)
		if (vertex.label.cex <= 0){ vertex.label.cex <- max(8*(nrow(sparseA))^(-0.6), 0.1) }
		if (vertex.size <= 0){ vertex.size <- max(75*(nrow(sparseA))^(-0.6), 1) }
		if (is.null(main)){ main <- "Partial correlation graph of VAR(1) model" }

		# reorder nodes for nicer plotting      
		shuffle <- maximum.cardinality.search(graph.adjacency(CIG, "undirected"))$alpham1
		CIG <- CIG[shuffle, shuffle]        
		nNames <- nNames[shuffle]        
		diag(CIG) <- 0
		gObj <- graph.adjacency(CIG, "undirected")

		# size of the edges       
		if (edge.width <= 0){ 
			edge.width <- 10 * (nrow(sparseA)^(-0.8)) 
		} else {
			edge.width <- edge.width 
		}
        
		# print(CIG)
		# layout=layout_(gObj, nicely())
		# actual plotting
		plot.igraph(gObj, layout=layout.circle, vertex.size=rep(vertex.size, nNodes), vertex.label.font=vertex.label.font, 
			vertex.color=rep(vertex.color.T0, nNodes), vertex.frame.color=vertex.frame.color,  edge.color="black",
			vertex.label=nNames, vertex.label.cex=vertex.label.cex, vertex.label.color=rep(vertex.label.color.T0, nNodes), 
			edge.width=edge.width, vertex.label.family="sans", main=main, ...)       
	}
    
	if (type=="contempPC"){
		# prune covariate without connections (according to sparseA and sparseP)
		if (prune){ 
			diag(sparseP) <- 0 
			idRemove <- intersect(which(apply(sparseP, 1, function(Z){ all(Z == 0) })), which(apply(sparseP, 2, function(Z){ all(Z == 0) })))
			if (length(idRemove) > 0){ 
				sparseA <- sparseA[-idRemove, -idRemove] 
				sparseP <- sparseP[-idRemove, -idRemove] 
				nNames <- nNames[-idRemove]
			}
		}
		nNodes <- nrow(sparseP)
		if (vertex.label.cex <= 0){ vertex.label.cex <- max(8*(nrow(sparseP))^(-0.6), 0.1) }
		if (vertex.size <= 0){ vertex.size <- max(75*(nrow(sparseP))^(-0.6), 1) }
		if (is.null(main)){ main <- "Contemp. cond. independence graph of VAR(1) model" }
 
		# get global CIG
		CIG <- CIGofVAR1(sparseA, sparseP, type="contemp")
		shuffle <- maximum.cardinality.search(graph.adjacency(CIG, "undirected"))$alpham1
		nNames <- nNames[shuffle]        
		CIG <- CIG[shuffle, shuffle]
		diag(CIG) <- 0
		gObj <- graph.adjacency(CIG, "undirected")

		# size of the edges       
		if (edge.width <= 0){ edge.width <- 10 * (nrow(sparseP)^(-0.8)) }
		if (edge.width > 0){ edge.width <- edge.width }
        
		# actual plotting
		plot.igraph(gObj, layout=layout.circle, vertex.size=rep(vertex.size, nNodes), vertex.label.font=vertex.label.font, 
                	vertex.color=rep(vertex.color.T0, nNodes), vertex.frame.color=vertex.frame.color,  edge.color="black",
                	vertex.label=nNames, vertex.label.cex=vertex.label.cex, vertex.label.color=rep(vertex.label.color.T0, nNodes), 
                	edge.width=edge.width, vertex.label.family="sans", main=main, ...)               
	}
}

