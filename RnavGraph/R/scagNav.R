##
## scagGraph functions
##
## These functions simplify the creation of navGraph sessions
## whose graphs (their nodes) follow certain scagnostic properties
###################################################################


## get scagnostic weights for scatterplots
scagEdgeWeights <- function(
    data,    ## data.frame or NG_data
    scags=c( # any subset list from 		
        "Clumpy", "NotClumpy",
        "Monotonic", "NotMonotonic",
        "Convex", "NotConvex",  
        "Stringy", "NotStringy",
        "Skinny", "NotSkinny",
        "Outlying","NotOutlying",
        "Sparse", "NotSparse",
        "Striated", "NotStriated",
        "Skewed", "NotSkewed"),
    combineFn = NULL  # function to combine across scagnostics
    ## Fraction of top values to use in (0,1]
) {
    ## Finds the combineFn of all scagnostics named
    ## NotSomething = 1 - Something
    ## Value is a weight vector for only those edges
    ## with combineFn of the named scagnostics
    ## having the highest topFrac of values.
    
    if(is(data,"NG_data")) {
        if(!(length(data@shortnames)==0)) {
            t.names <- shortnames(data)
            data <- data@data
            names(data) <- t.names
        }else {
            data <- data@data
        }
    }
    
    
    varNames <- names(data)
    
    scag <- scagnostics(data)  ## returns a matrix
    
    
    isNotSomething <- grepl("^not", tolower(scags))
    ## get the scagnostic strings
    scagNames <- gsub("^not","",tolower(scags))
    
    if(!all(scagNames %in%
            tolower(rownames(scagnostics(data.frame(
                a = c(1,0,1), b = c(2,3,0), c = c(1,9,3))))))) {
        stop('[scagEdgeWeights] argument scags contains invalid names.')
    }
    
    ## create empty from to edges matrix
    nNodes <- dim(data)[2]
    fromN <- rep(1:(nNodes-1), (nNodes-1):1)
    nE <- nNodes*(nNodes-1)/2
    toN <- rep(-99, nE)
    j <- 1
    for(i in seq(2,nNodes)){
        k <- nNodes-i+j
        toN[j:k] <- i:nNodes
        j<-k+1
    }
    
    ftE <- cbind(fromN,toN)
    colnames(ftE) <- c("from", "to")	
    
    
    ## make a matrix with the weights for each scags element
    iScags <- match(scagNames,tolower(rownames(scag)))
    
    scagColNames <- tolower(colnames(scag))
    
    weightsScag <- matrix(rep(-99, nE*length(scags)), nrow = nE)
    for(i in 1:nE) {
        name1 <- tolower(paste(varNames[ftE[i,1]],"*",varNames[ftE[i,2]]))
        name2 <- tolower(paste(varNames[ftE[i,2]],"*",varNames[ftE[i,1]]))
        
        i1 <- match(name1,scagColNames)
        i2 <- match(name2,scagColNames)
        
        if(is.na(i2)) {
            if(length(i1)!=1 && !is.na(i2)){
                stop("[scagEdgeWeights] non unique data names")
            }else {
                ii <- i1
            }
        }else {
            if(length(i2)!=1 && is.na(i1)){
                stop("[scagEdgeWeights] non unique data names")
            }else {
                ii <- i2
            }
        }
        
        ## incorporate NotSomething
        weightsScag[i,] <- isNotSomething + scag[iScags,ii] * (-(isNotSomething*2-1))		
    }
    colnames(weightsScag) <- scags
    
    
    if(is.null(combineFn)) {
        return(list(fromToEdgeMatrix = cbind(ftE,weightsScag), nodeNames = varNames))
        
    }else {		
        FUN <- match.fun(combineFn)
        ftEmat <- cbind(ftE,apply(weightsScag,1,FUN))
		colnames(ftEmat) <- c("from","to", "combined weights")
        
        return(list(fromToEdgeMatrix = ftEmat, nodeNames = varNames))
    }
    
    
}


## create variable graphs according to scagnostic measures
scagGraph <- function(edgeWeights, topFrac = 0.2){
    
    if(dim(edgeWeights$fromToEdgeMatrix)[2]==3){
        weights <- edgeWeights$fromToEdgeMatrix[,3] 
        if(topFrac != 1){
            ii <- weights>quantile(weights,1-topFrac)
        }else {
            ii <- rep(TRUE,length(weights))
        }
        if(all(ii == FALSE)) {
            return(newgraph(nodeNames = edgeWeights$nodeNames,
                            mat = matrix(rep(0,length(edgeWeights$nodeNames)^2),
                                ncol = length(edgeWeights$nodeNames)),isAdjacency=TRUE))
        }else {			
            return(newgraph(nodeNames = edgeWeights$nodeNames,
                            mat = edgeWeights$fromToEdgeMatrix[ii,c(1,2)],
                            weights = weights[ii]))
        }
    }else {
        ftE <- edgeWeights$fromToEdgeMatrix[,c(1,2)]
        scagWeights <- edgeWeights$fromToEdgeMatrix[,-c(1,2)]
        
        nodeNames <- edgeWeights$nodeNames
        
        nCol <- dim(scagWeights)[2]
        graphList <- vector("list" , length = nCol)
        
        names(graphList) <- colnames(scagWeights)
        
        for(i in 1:nCol){
            weights <- scagWeights[,i]
            if(topFrac != 1){
                ii <- weights>quantile(weights,1-topFrac)
            }else {
                ii <- rep(TRUE,length(weights))
            }
            if(all(ii == FALSE)) {
                graphList[[i]] <- newgraph(nodeNames = edgeWeights$nodeNames,
                                           mat = matrix(rep(0,length(edgeWeights$nodeNames)^2),
                                               ncol = length(edgeWeights$nodeNames)))
            }else {
                graphList[[i]] <- newgraph(nodeNames = edgeWeights$nodeNames,
                                           mat = cbind(ftE[ii,],scagWeights[ii,i]), weights = weights[ii])				
            }
        }
        return(graphList)
    }
}





##
##  Putting it all together with navGraph
##
scagNav <- function (data,   # ngdata object
                     scags=c(	# any subset list from 
                         "Clumpy", "NotClumpy",
                         "Monotonic", "NotMonotonic",
                         "Convex", "NotConvex",  
                         "Stringy", "NotStringy",
                         "Skinny", "NotSkinny",
                         "Outlying","NotOutlying",
                         "Sparse", "NotSparse",
                         "Striated", "NotStriated",
                         "Skewed", "NotSkewed"),
                     topFrac=0.2,
                     ## Fraction of top values to use in (0,1],
                     combineFn=NULL,
                     ## NULL means all of the scags are done
                     ## and laid out in the graphs menu 
                     ## Otherwise, combineFn should be any
                     ## *NAMED* function on a matrix which
                     ## returns a single numerical value, 
                     ## e.g. max, sum, or one of your own 
                     ## making.
                     settings = NULL,		# Settings for navgraph - a list of name value settings.
                     glyphs = NULL,
                     images = NULL,
                     sep = ":",
                     layout="circle"  
                     ## layout is identical to that for ng_graph
                     ## "circle", 
                     ## "random", 
                     ## "kamadaKawaiSpring",
                     ## "fruchtermanReingold"
                     ) {
    edgeWeights <- scagEdgeWeights(data = data, scags = scags, combineFn = combineFn)
    
    graphList <- scagGraph(edgeWeights, topFrac = topFrac)
    
    if(length(graphList) == 1){
        graphList <- list(graphList)
    }
    
    ng.graphList <- vector("list", length = 2 * length(graphList))
    ng.vizList <- vector("list", length = 2 * length(graphList))
    
    for(i in seq(1,2*length(graphList),by = 2)) {

        j <- ceiling(i/2)
        LG <- linegraph(graphList[[j]], sep = sep)
        ng.graphList[[i]] <- ng_graph(name = paste("3d", names(graphList)[j]),
                                      graph = LG, sep = sep,layout = layout)
        ng.graphList[[i+1]] <- ng_graph(name = paste("4d", names(graphList)[j]),
                                        graph = complement(LG), sep = sep,
                                        layout = layout)
        ng.vizList[[i]] <- ng_2d(data, ng.graphList[[i]],
                                 glyphs = glyphs, images = images)
        ng.vizList[[i+1]] <-  ng_2d(data, ng.graphList[[i+1]],
                                    glyphs = glyphs, images = images)
    }
    
    ##
    ## Set up the navGraph
    ng <- navGraph(data, ng.graphList, ng.vizList, settings = settings)
    
    ## return ng handler
    return(ng)
}



