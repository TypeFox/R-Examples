setMethod("initialize", signature(.Object="RSeed"),
	function(.Object, model, connectedComponentCutOff=10, currencyMetabolites){
		packs <- c(sybil=require(sybil), graph=require(graph), RBGL=require(RBGL))
		
		.Object@env <- new.env()
		
		if(!all(packs)){
			print(packs)
			stop("missing packages")
		}
		
		if(class(model) != "modelorg"){
			stop("need modelorg")
		}
		
		if(!missing(currencyMetabolites)){
			currencyMetabolites(.Object) <- currencyMetabolites
		}
		
		connectedComponentCutOff(.Object) <- connectedComponentCutOff
		model_original(.Object) <- model
		return(.Object)
	}
)


RSeed <- function(model=NULL, connectedComponentCutOff=10){
	b <- new("RSeed", model, connectedComponentCutOff)
	
	return(b)
}
setGeneric("makeExperiment",
		function (rs) standardGeneric("makeExperiment")
)
setMethod("makeExperiment", signature(rs="RSeed"),
	function(rs){
		buildGraph(rs)
		sourceCompounds(rs)
		return(rs)
	}
)
	
setGeneric("buildGraph",
		function (rs) standardGeneric("buildGraph")
)
setMethod("buildGraph", signature(rs="RSeed"),
	function(rs){

		met_rm <- NULL
		rmReact <- NULL
	
		cat("building graph\n")
		
		model <- model_original(rs)
		
		mygraphAM <- new("graphAM", edgemode="directed", adjMat=buildAdjacency(model, currencyMetabolites(rs)))
	
		#strongComp is faster on graphNEL
		mygraph <- as(mygraphAM, Class="graphNEL")
		rm(mygraphAM)
	
		cat("preprocessing: checking if graph is connected\n")
	
		cc <- connectedComp(mygraph)
	
		whichcc <- which(sapply(cc, length) < connectedComponentCutOff(rs))
	
		if(length(cc) > 1){
		cat("graph is not connected, trying to fix this\n")
	
			if( length(cc) - length(whichcc) > 2){
				stop("your graph has not connected components, which are bigger than \'connectedComponentCutOff\' nodes")
			}
		
			cat(paste("removing parts smaller than ", connectedComponentCutOff(rs), " nodes from model\n", sep="", collapse=""))
		
			met_rm <- character(0)
		
			##TODO look for elegant way to solve.
			for(i in cc[whichcc]){
				met_rm <- c(met_rm, i)
			}
			
	#		pos <- sapply(met_rm,	function(x){
	#									grep(pattern=x, x=met_id(model))
	#								}
	#		)

			pos <- as.integer(sub("n", "", met_rm))


#			dont need this			
#			if(!is.null(currencyMetabolites(rs))){
#				#pos <- c(pos,which(met_id(model) %in% currencyMetabolites(rs)))
#			}
			
			
			cat("following metabolites were removed\n")
			met_rm <- met_id(model)[pos]
			print(met_rm)
		
			#TODO entfernt wohl nicht ganz sauber
			met_name(model) <- met_name(model)[-pos]
			met_id(model) <- met_id(model)[-pos]
			met_comp(model) <- met_comp(model)[-pos]
			met_single(model) <- met_single(model)[-pos]
			rhs(model) <- rhs(model)[-pos]

			met_num(model) <- met_num(model) - length(pos)
		
			S(model) <- S(model)[-pos, ]
		
		
			cat("checking for reactions without metabolites and remove them\n")
		
			pos <- which(apply((S(model)==0), 2, all))
		
			cat("following reactions were removed:\n")
		
			print(react_id(model)[pos])
			rmReact <- react_id(model)[pos]
			
			model <- rmReact(model, pos)
			cat("not connected parts removed, building new graph\n")
		
			model_edited(rs) <- model
		
			mygraphAM <- new("graphAM", edgemode="directed", adjMat=buildAdjacency(model, currencyMetabolites(rs)))
			#strongComp is faster on graphNEL
			mygraph <- as(mygraphAM, Class="graphNEL")
			rm(mygraphAM)
		}
		else{
			cat("graph is connected\n")
		}
	
		graph_network(rs) <- mygraph
		model_changes(rs) <- list(rmMet=met_rm, rmReact=rmReact)
	
		scc(rs) <- strongComp(mygraph)
})


setGeneric("sourceCompounds",
		function (rs) standardGeneric("sourceCompounds")
)


setMethod("sourceCompounds", signature(rs = "RSeed"),
	function(rs){
		validObject(rs)
		
#		if(is.null(graph)){
#			stop("need graph object")
#		}
#		if(is.null(sc)){
#			stop("need strong components as a list")
#		}
#		if(class(graph)!="graphNEL"){
#			stop("graph has to be a graph")
#		}
#		if(!is.list(sc)){
#			stop("sc has to be a list")
#		}
		
		g <- graph_network(rs)
		scc <- scc(rs)
		
		scc_sizes <- sapply(scc, length)

		n <- sapply(seq(along=scc(rs)), 
			function(x){
				if(scc_sizes[x] > 1){
					return(paste("c", as.character(x), sep="", collapse=""))
				}
				else{
					return(scc[[x]][1])
				}
			}
		)

		names(scc(rs)) <- n
		scc <- scc(rs)
		names(scc_sizes) <- n
		scc_sizes(rs) <- scc_sizes

		lookup <- unlist(lapply(names(scc), 
			function(x){ a <- rep(x, scc_sizes[x])
			names(a) <- scc[[x]]
			return(a)
		}))
		
		is_seed <- rep(TRUE, length(scc))
		names(is_seed) <- names(scc)
		
		
		## build overlay-graph
		am <- matrix(0, ncol=length(scc), nrow=length(scc), dimnames=list(names(scc), names(scc)))
		for(x in nodes(g)){
			lx <- lookup[x]
			for(ly in lookup[edges(g, x)[[1]]]){
				#ly <- lookup[y]
				if(lx != ly){
					
					am[lx, ly] <- 1
				
#						if(!isAdjacent(g_scc, from=lx, to=ly)){
#							g_scc <- addEdge(graph=g_scc, from=lx, to=ly)
#						}
				}
			}
			
		}
		
		g_scc <- new("graphAM", adjMat=am, edgemode="directed")
		
		graph_scc(rs) <- g_scc
		list_sc(rs) <- names(which(degree(g_scc)$inDegree == 0))
		
#		if(FALSE){
#		t2 <- system.time({
#		for(x in scc(rs)){ 
#			# take every node from the scc
#			for(i in x){
#				# take every edge of this node
#				
#				e <-  edges(g, i)[[1]]
#				
#				drin <- e %in% x
#				
#				is_seed[lookup[e[!drin]]] <- FALSE
#				
#				for(j in edges(g, i)[[1]]){
#					# is this edge directed to a node in this scc (i)?
#					if(length(j)==0){
#						next;
#						cat("blub\n")
#					}
#					if(!any(j == x)){
#				
#						# scc, which contains j is not a seed scc!
#						is_seed[lookup[j]] <- FALSE
#					}
#				}
#			}
#		}
#
#		list_sc(rs) <- names(which(is_seed))
#		})
#		}
#		return(list(t1, t2))
		
		
		
		#list_sc
		#scc_sizes
	}
)


setGeneric("confidenceLevel",
		function (rs) standardGeneric("confidenceLevel")
)

setMethod("confidenceLevel", signature(rs = "RSeed"),
	function(rs){
		level <- 1/scc_sizes(rs)[list_sc(rs)]
		match <- level[level >= 0.2]
		return(list(confidenceLevel=level, match_threshold=match))
	}
)




setGeneric("findCurrencyMetabolites",
		function (object, ...) standardGeneric("findCurrencyMetabolites")
)

setMethod("findCurrencyMetabolites", signature(object = "modelorg"),
	function(object, cutoff=20){
		
		cm <- sybil::met_id(object)[which(Matrix::rowSums(x=(sybil::S(object)!= 0) ) > cutoff)]
		return(cm)
	}
)


setMethod("plot", signature(x="RSeed", y="missing"),
	function(x, ...){
		if(!isTRUE(require("Rgraphviz"))){
			stop("need Rgraphviz")
		}
		
		rs <- x
		
		if(any(is.null(scc(rs)), is.null(list_sc(rs)))){
			stop("no results to plot, first run \"makeExperiment\"")
		}
		
		## this part should not be neccessary any more
		if(is.null(graph_scc(rs))){
			cat("aggregating nodes, plz be patient\n")
			g_scc <- graph_network(rs)

			for(i in which(scc_sizes(rs)>1)){
				g_scc <- combineNodes(scc(rs)[[i]], g_scc, names(scc(rs)[i]))
			}
			
			graph_scc(rs) <- g_scc
		}
		
		layout(matrix(data = 1:3, nrow = 1, ncol = 3), 
				widths = c(8, 1, 1), heights = c(1))
		
		len <- sapply(scc(rs), length)
		
		names(len) <- names(scc(rs))
		
		len <- cbind(len, 0)
		len[list_sc(rs), 2] <- 1
		
		isSeed <- len[,2] == 1
		
		uSeed <- unique(len[isSeed,1])
		unSeed <- unique(len[!isSeed,1])
		
		seedcolors <- rainbow(length(uSeed), start=0.3, end=0.7)
		names(seedcolors) <- uSeed
		
		nseedcolors <- rainbow(length(unSeed), start=0.8, end=0.2)
		names(nseedcolors) <- unSeed
		
		
		colors <- apply(len, 1,	function(x){
									if(x[2]==1){
										return(seedcolors[as.character(x[1]) ])
									}
									else{
										return(nseedcolors[as.character(x[1])])
									}
		
								}
		)
		
		vh <- rep(20, length(nodes(graph_scc(rs))))
		names(vh) <- names(colors)
		
		nattrs <- list(fillcolor=colors, height=vh)
		gattrs <- list(node = list(shape = "box", fixedsize = FALSE))
		
		dn <- function(node, ur, attrs, radConv) {
			nodeCenter <- getNodeCenter(node)
			nodeX <- getX(nodeCenter)
			nodeY <- getY(nodeCenter)
			lw <- getNodeLW(node)
			rw <- getNodeRW(node)
			rad <- (lw + rw)/2
			height <- getNodeHeight(node)
			fg <- color(node)
			style <- style(node)
			if (fg == "") 
				fg <- "black"
			bg <- fillcolor(node)
			if (bg == "") {
				if (style == "filled") 
				    bg <- "grey"
				else bg <- "transparent"
			}
			rect(nodeX - lw, nodeY - (height/2), nodeX + rw, nodeY + (height/2), col = bg, border = fg)
		
			#drawTxtLabel(txtLabel(node), xLoc = nodeX, yLoc = nodeY)
	
	
			text(name(node), x=getX(nodeCenter), y=getY(nodeCenter), srt=90, cex=1)
		}


		plot(agopen(graph_scc(rs), ..., nodeAttrs=nattrs, attrs=gattrs, name=""), drawNode=dn)
		title(mod_name(model_original(rs)))
		
#		pie(1:10, col=rainbow(10))
		
		par(mar = c(5,2,5,2))
		
		v1 <- 1:length(uSeed)
		m <- matrix(NA, nrow=9, ncol=length(uSeed))
		m[5,] <- v1
		
		image(x=1:9, y=v1, z=m,
				 col=seedcolors[order(uSeed)], ylab="", xlab="", xaxs="i", yaxs="i", xaxt='n', yaxt='n', bty='n')
		axis(side=2, at=v1, labels=sort(uSeed))
		title("Seeds")
		
		
		v2 <- 1:length(unSeed)
		m <- matrix(NA, nrow=9, ncol=length(unSeed))
		m[5,] <- v2
		
		image(x=1:9, y=v2, z=m,
				 col=nseedcolors[order(unSeed)], ylab="", xlab="", xaxt='n', yaxt='n', bty='n')
		axis(side=2, at=v2, labels=sort(unSeed))
		title("notSeeds")
		
		
		
		layout(1)
	}
)


setGeneric("getSourceMetabolites",
		function (rs) standardGeneric("getSourceMetabolites")
)

setMethod("getSourceMetabolites", signature(rs = "RSeed"),
	function(rs){
		sourceCompounds <- scc(rs)[list_sc(rs)]
		sourceCompounds <- lapply(sourceCompounds,
			function(x){
				met_id(model_used(rs))[as.integer(gsub("n", "", x))]
			}
		)
		return(sourceCompounds)
	}
)





buildAdjacency <- function(model="modelorg", currencyMetabolites=NULL){

	if(class(model) != "modelorg"){
		stop("need modelorg")
	}
	
	graphpackage <- require("graph")
	
	if(!isTRUE(graphpackage)){
		stop("need graphpackage for this job")
	}
	
	adj <- matrix(0, nrow=met_num(model), ncol=met_num(model))#, dimnames=list(met_id(model), met_id(model)))
	
	apply(S(model), 2,
		function(x){
			source <- which(x > 0)
			#dest <- which(x < 0)
			return(source)
#			if(length(source) + length(dest) > 4){
#				stop("reaction with more than 4 metabolites")
#			}
			
		}
	)
	
	for(i in seq(along=react_id(model))){
		source <- which(S(model)[,i] < 0)
		dest <- which(S(model)[,i] > 0)
		
		adj[source, dest] <- 1
		
		if(react_rev(model)[i]){
			adj[dest, source] <- 1
		}
	
	}
	
	if(!is.null(currencyMetabolites)){
		met <- which(met_id(model) %in% currencyMetabolites)
		adj[, met] <- 0
		adj[met,] <- 0
		
	}
	
	return(adj)
}
