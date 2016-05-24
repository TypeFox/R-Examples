########### optCluster Function ###########

## Determine optimal clustering method and number of clusters
optCluster <- function(obj, nClust, clMethods = "all", countData = FALSE, validation = c("internal", "stability"), hierMethod = "average", 
		 annotation = NULL, clVerbose = FALSE, rankMethod = "CE", distance = "Spearman", importance = NULL, rankVerbose= FALSE, ...) {

	is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol     
    isCountData <- all(is.wholenumber(obj))
    
    if(isTRUE(isCountData) && identical(countData, FALSE)){
    	warning("obj is count data but countData is set to FALSE")
    } else if(identical(isCountData, FALSE) && isTRUE(countData)){
    	warning("obj is not count data but countData is set to TRUE")
    }
     	
	if("all" %in% clMethods) {
		if(isTRUE(countData)){
			clMethods <- c("em.nbinom", "da.nbinom", "sa.nbinom", "em.poisson", "da.poisson", "sa.poisson")
		} else {
			if(hierMethod == "ward"){
			clMethods <- c("agnes", "clara", "diana", "fanny", "hierarchical", "kmeans", "model", "pam", "som", "sota")
			} else{
			clMethods <- c("clara", "diana", "fanny", "hierarchical", "kmeans", "model", "pam", "som", "sota")
			}
		}
    }
  	
  	if("all" %in% validation) {
	validation <- c("internal", "stability", "biological")
    }  	
	
	nClustOut <- nClust
	switch(class(obj), matrix = mat <- obj, ExpressionSet = mat <- Biobase::exprs(obj), 
        data.frame = {
            if (any(!sapply(obj, class) %in% c("numeric", "integer"))) stop("data frame 'obj' contains non-numeric data")
            mat <- as.matrix(obj)
        }, stop("argument 'obj' must be a matrix, data.frame, or ExpressionSet object"))	
		
	addArgs <- list(...)
	## 'verbose' and 'method' are used in more than one function
	if(exists(c("verbose"), where = addArgs)){
		stop("must specify 'verbose' as 'clVerbose' or 'rankVerbose'")
	}
	if(exists(c("method"), where = addArgs)){
		stop("must specify 'method' as 'hierMethod', or 'rankMethod'")
	}
	mbCountNames <- c("em.nbinom", "da.nbinom", "sa.nbinom", "em.poisson", "da.poisson", "sa.poisson")
	mbCountMethods <- clMethods[clMethods %in% mbCountNames]	
	
	if(length(mbCountMethods) > 0){
		if(!requireNamespace("MBCluster.Seq")) {
      		stop("package 'MBCluster.Seq' required for clustering count data")
    	}
   	
     	if(isTRUE(isCountData)){	
		mbCount <- TRUE
		clMethods <- clMethods[-which(clMethods %in% mbCountNames)]
		} else {
			stop("obj is not count data and or more of the model based clustering algorithms for count data have been selected")
		}
		## remove any additional arguments related to clustering algorithms for count data
		if(exists(c("Treatment"), where = addArgs)){
			uniqueTreatment <- addArgs$Treatment
			addArgs <- addArgs[-which(names(addArgs) == "Treatment")]
		} else if(is.null(colnames(obj))) {
			stop("column names or 'Treatment' argument must be used for clustering count data")
		} else {
			uniqueTreatment <- colnames(obj)
			if(length(uniqueTreatment) == length(unique(uniqueTreatment))){
				## find unique treatments
				uniqueTreatment <- substr(uniqueTreatment,1,nchar(uniqueTreatment)-1)		
			}
		}
		if(exists(c("Normalizer"), where = addArgs)){
			Normalizer <- addArgs$Normalizer
			addArgs <- addArgs[-which(names(addArgs) == "Normalizer")]
		} else{
			Normalizer <- NULL
			}
		if(exists(c("iter.max"), where = addArgs)){
			iter.max <- addArgs$iter.max
			addArgs <- addArgs[-which(names(addArgs) == "iter.max")]
		} else{
			iter.max <- 30
			}
		if(exists(c("TMP"), where = addArgs)){
			TMP <- addArgs$TMP
			addArgs <- addArgs[-which(names(addArgs) == "TMP")]
		} else{
			TMP <- NULL
		}					
	} else {
		mbCount <- FALSE
	}
	
	## Sort arguments to clValid or RankAggreg	
	RankAggreg.names <- c(names(formals(RankAggreg)),"p")
	RankAggreg.args <- addArgs[names(addArgs) %in% RankAggreg.names]
	if(length(RankAggreg.args) > 0){
		clValid.args <- addArgs[-which(names(addArgs) %in% RankAggreg.names)]
	} else{
		clValid.args <- addArgs	
	}
	
	if(exists(c("maxitems"), where = clValid.args)){
		maxitems <- clValid.args$maxitems
		clValid.args <- clValid.args[-which(names(clValid.args) == "maxitems")]
	}else if(nrow(obj) > 600){
		maxitems <- nrow(obj)
	}else {
		maxitems <- 600
	}
	
	## Obtain validation measures for clValid clustering algorithms
	if(length(clMethods) > 0){
		clVObj <- do.call('clValid',c(list(obj, nClust, clMethods, validation, method = hierMethod, annotation = annotation, verbose = clVerbose,
				maxitems = maxitems), clValid.args))
		if(any(is.na(clVObj@measures))){
			stop("NA validation measures returned for at least one clustering algorithm \n  Rank aggregation cannot be performed")
		}
	  	clValResults <- clVObj@clusterObjs
	  	measNames <- measNames(clVObj)
  		clAlgs <- clusterMethods(clVObj)
  		metric <- clVObj@metric
  		neighbSize <- clVObj@neighbSize
  		GOcategory <- clVObj@GOcategory
  		goTermFreq <- clVObj@goTermFreq
  	} else {
  		## Obtain arguments for calculating only model based clustering for count data
  		measNames <- c(if("stability"%in%validation) c("APN","AD","ADM","FOM"),
                 if("internal"%in%validation) c("Connectivity","Dunn","Silhouette"),
                 if("biological"%in%validation) c("BHI","BSI"))
        clValResults <- list()
        clAlgs <- vector()
        if(exists(c("metric"), where = clValid.args)){
			metric <- clValid.args$metric
		} else{
			metric <- "euclidean"
		}
		if(exists(c("neighbSize"), where = clValid.args)){
			neighbSize <- clValid.args$neighbSize
		} else{
			neighbSize <- 10
		}	
		
		if(exists(c("GOcategory"), where = clValid.args)){
				GOcategory <- clValid.args$GOcategory
			} else {
				GOcategory <- "all"
		}
		if(exists(c("goTermFreq"), where = clValid.args)){
				goTermFreq <- clValid.args$goTermFreq
			} else {
				goTermFreq <- 0.05
		}
		if(exists(c("dropEvidence"), where = clValid.args)){
				dropEvidence <- clValid.args$dropEvidence
		} else {
				dropEvidence <- NULL
		}		
	}       
  
  		
	## Obtain validation and clustering results for model based clustering for count data
    if(isTRUE(mbCount)){
    	valMBCount <- mbCountVal(mat, Normalizer, uniqueTreatment, validation, measNames, nClust, metric, neighbSize, clVerbose, 
 					mbCountMethods, iter.max, TMP, annotation, GOcategory, goTermFreq, dropEvidence)
		MBCountArray <- valMBCount$measures
		MBCountCluster <- valMBCount$clusterObj	
	}		
	## Merge mbSeq results with clValid results
  	nClust <- as.character(nClust)
  	if(isTRUE(mbCount)){
  		clAlgs <- c(clAlgs, mbCountMethods)
  		meas <- array(dim = c(length(measNames), length(nClust), length(clAlgs)))
  		dimnames(meas) <- list(measNames,nClust,clAlgs)
  		j <- 1
  		for(i in 1:dim(meas)[3]){
  			if(i <= (dim(meas)[3]-length(mbCountMethods))){
  				meas[,,i] <- clVObj@measures[,,i]	
  			} else {
  				meas[,,i] <- MBCountArray[,,j]
  				j <- j+1	
  			}	
  		} 		
  		clusterResults <- c(clValResults, MBCountCluster)
   	} else{
  		meas <- clVObj@measures
  		clusterResults <- clValResults
  	}
      
	## Create new 'clValid' class object
	clVal <- new("clValid", clusterObjs = clusterResults, measures = meas, 
        measNames = measNames, clMethods = clAlgs, labels = rownames(obj), 
        nClust = nClustOut, validation = validation, metric = metric, 
        method = hierMethod, neighbSize = neighbSize, GOcategory = GOcategory, 
        goTermFreq = goTermFreq, annotation = annotation, call = match.call())
    
    ## Obtain clustering ranks and weights    
 	clusterOrder <- getRanksWeights(clVal)
 	if("biological" %in% validation){
	BHIrowNum <- which(rownames(clusterOrder$ranks) == "BHI")
	clusterOrder$ranks[BHIrowNum,] <- rev(clusterOrder$ranks[BHIrowNum,])
	clusterOrder$weights[BHIrowNum,] <- rev(clusterOrder$weights[BHIrowNum,])	
	BSIrowNum <- which(rownames(clusterOrder$ranks) == "BSI")
	clusterOrder$ranks[BSIrowNum,] <- rev(clusterOrder$ranks[BSIrowNum,])
	clusterOrder$weights[BSIrowNum,] <- rev(clusterOrder$weights[BSIrowNum,])
	}
		
	if(exists(c("k"), where = RankAggreg.args)){
		nList <- RankAggreg.args$k
		RankAggreg.args <- RankAggreg.args[-which(names(RankAggreg.args) == "k")]
	} else {
		nList <- ncol(clusterOrder$ranks)
	}
	if(exists(c("weights"), where = RankAggreg.args)){
		clusterOrder$weights <- RankAggreg.args$weights
		RankAggreg.args <- RankAggreg.args[-which(names(RankAggreg.args) == "weights")]
	}
	## Weights for validation measure lists
 	if(is.null(importance)){
 	  importance <- rep(1, nrow(clusterOrder$ranks))
 	}
 	
	## Perform weighted rank aggregation		
	optimal.list <- do.call('RankAggreg', c(list(x = clusterOrder$ranks, k = nList, weights = clusterOrder$weights, 
	method = rankMethod, distance = distance, importance = importance, verbose = rankVerbose), RankAggreg.args))
        	
	## Create 'optCluster' class object
	new("optCluster", inputData = mat, clVal = clVal, ranksWeights = clusterOrder, rankAgg = optimal.list)	
}

########## repRankAggreg Function ##########

## Repeat weighted rank aggregation on "optCluster" object
repRankAggreg <- function(optObj, rankMethod = "same", distance = "same", importance = NULL, rankVerbose = FALSE, ... ){
	
	## Obtain RankAggreg arguments
	RankAggreg.args <- list(...)
	clusterOrder <- list(ranks = methodRanks(optObj), weights = scoreRanks(optObj))
	rankAgg <- getRankAggreg(optObj)
	method <- match.arg(rankMethod, c("same", "CE", "GA"))
	distance <- match.arg(distance, c("same", "Spearman", "Kendall"))
	
	if(method == "same") {
		method <- rankAgg$method	
	} 
	if(distance == "same") {
		distance <- rankAgg$distance	
	} 
	if(exists(c("method"), where = RankAggreg.args)){
		message(" The argument 'method' has been changed to 'rankMethod' ")
		method <- RankAggreg.args$method
		RankAggreg.args <- RankAggreg.args[-which(names(RankAggreg.args) == "method")]	
	} 		
	if(exists(c("verbose"), where = RankAggreg.args)){
		message(" The argument 'verbose' has been changed to 'rankVerbose' ")
		rankVerbose <- RankAggreg.args$verbose
		RankAggreg.args <- RankAggreg.args[-which(names(RankAggreg.args) == "verbose")]	
	} 		
	if(exists(c("k"), where = RankAggreg.args)){
		nList <- RankAggreg.args$k
		RankAggreg.args <- RankAggreg.args[-which(names(RankAggreg.args) == "k")]
	} else {
		nList <- ncol(clusterOrder$ranks)
	}
	if(exists(c("weights"), where = RankAggreg.args)){
		clusterOrder$weights <- RankAggreg.args$weights
		RankAggreg.args <- RankAggreg.args[-which(names(RankAggreg.args) == "weights")]
	}
	
	## Weights for validation measure lists
	if(is.null(importance)){
	  importance <- rep(1, nrow(clusterOrder$ranks))
	}
	
	## Perform weighted rank aggregation		
	optimal.list <- do.call('RankAggreg', c(list(x = clusterOrder$ranks, k = nList, weights = clusterOrder$weights, 
	method = method, distance = distance, importance = importance, verbose = rankVerbose), RankAggreg.args))
	
	## Create new 'optCluster' class object
		new("optCluster", inputData = getDataset(optObj), clVal = getClValid(optObj), ranksWeights = clusterOrder, rankAgg = optimal.list)
}

########## Plot Functions for optCluster Object ##########

## Validation Measure Plots
valPlot <- function(x, measures=measureNames(x), legend=TRUE, legendLoc="topright", main=NULL,
                   pch=NULL, type="b", ask=prod(par("mfcol")) < length(measures) && dev.interactive(), ...) {
          	clValObject <- getClValid(x)
          	plot(clValObject, measures = measures, legend = legend, legendLoc = legendLoc, main = main,
          		pch = pch, type = type, ask = ask, ...)
}  
          

## Rank Aggregation Plots
aggregPlot <- function(x, show.average = TRUE, show.legend = TRUE, colR = "red",...) {
			raggr <- getRankAggreg(x)
			plot(raggr, show.average = show.average, show.legend = show.legend, colR = colR, ... )
} 
          
## Hierarchical Heat Map
optHeatmap <- function(x, dendroClusters = TRUE, barClusters = FALSE, clusterColors = "rainbow",
					mapColors = colorRampPalette(c("green", "black", "red"))(256), Colv = FALSE,
					dendrogram = "row", density.info = "none", ...){
			if(!requireNamespace("gplots")) {
      			stop("package 'gplots' required for plotting heat map")
    		}
    		## Obtain arguments for heatmap.2
    		plotArgs <- list(...)
    		if(exists("trace", where = plotArgs)){
    			trace <- plotArgs$trace
    			plotArgs <- plotArgs[-which(names(plotArgs) == "trace")]
    		} else {
    			trace <- "none"    		
    		} 
 
    		if(exists("scale", where = plotArgs)){
    		  scale <- plotArgs$scale
    		  plotArgs <- plotArgs[-which(names(plotArgs) == "scale")]
    		} else {
    		  scale <- "none"    		
    		} 
    		   		
    		## Obtain optimal clustering method and number of clusters 		    		
			topAlg <- topMethod(x)
			k <-as.numeric(gsub("\\D", "", topAlg))
			algName <- as.character(gsub("\\d", "", topAlg))
			algName <- gsub("-", "", algName)
			 
			## Obtain clustering results
			hierAlgs <- c("hierarchical", "agnes", "diana")
			if(algName %in% hierAlgs){
				res <- clusterResults(x, algName)
				den <- as.dendrogram(res)
				assignment <- cutree(res,k)			
			} else{
				stop("Optimal clustering algorithm is not hierarchical, agnes, or diana")
			}
			
			## Choose colors for clusters
			if(length(clusterColors) == 1 && clusterColors == "rainbow"){
				colors <- rainbow(k)
			}else{
				colors <- clusterColors[1:k]
			}			
			
			## Color dendrogram clusters
			if(isTRUE(dendroClusters)){
				orderh <- sort(res$height, decreasing = TRUE)
				minh <- orderh[k]
				maxh <- orderh[k-1]
				cuth <- mean(minh,maxh)
				denCut <- cut(den,cuth)
				hlist <- sort(orderh[1:k])
				for(j in 1:k){
					if(!is.element(attr(denCut$lower[[j]], "height"),hlist) && attr(denCut$lower[[j]], "height") > 0){
						hlist[length(hlist)+1] <- attr(denCut$lower[[j]], "height")					
					} 
				}
				colbranches <- function(n) {
					a <- attributes(n)
					if(!is.element(a$height, hlist)){					
						colN <- which(rownames(getDataset(x)) == labels(n)[1])
						attr(n, "edgePar") <- c(a$edgePar, list(col=colors[assignment[colN]], lwd=2))
					}
					n
				}
				Rowv <- dendrapply(den, colbranches)
			}else{
				if(exists(c("Rowv"), where = plotArgs)){
				Rowv <- plotArgs$Rowv
				plotArgs <- plotArgs[-which(names(plotArgs) == "Rowv")]
				}else {
				Rowv <- den					
				}
			}
			
			## Color sidebar clusters and plot heat map
			if(isTRUE(barClusters)){
			RowSideColors <- colors[assignment]
			do.call(heatmap.2, c(list(getDataset(x), Rowv = Rowv, Colv = Colv, dendrogram = dendrogram, scale = scale,
			 trace = trace, density.info = density.info, col = mapColors, RowSideColors = RowSideColors), plotArgs))
			}
			else{
			do.call(heatmap.2, c(list(getDataset(x), Rowv = Rowv, Colv = Colv, dendrogram = dendrogram, scale = scale,
			 trace = trace, density.info = density.info, col = mapColors), plotArgs))			
			}				
}
