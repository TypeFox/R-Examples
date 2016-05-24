setClass(
  Class = "NG_graph",
  representation = representation(
    name = "character",
    graph = "graph",
    xNodes = "numeric",
    yNodes = "numeric",
    xLabels = "numeric",
    yLabels = "numeric",
    bbox = "numeric",
    layout = "character",
    border = "numeric",
    sep = 'character'
    ),
  validity = function(object){
    goThrough <- TRUE
    
    if(is(object@graph,'graph') == FALSE){
      cat('[NG_graph:validation] graph is not from class graph.\n')
      goThrough <- FALSE
    }
    
    if( !all(numNodes(object@graph) == c(length(object@xNodes),length(object@yNodes),
                       length(object@xLabels), length(object@yLabels), length(object@xNodes)))){
      cat('[NG_graph:validation] dimension of a node related slot is wrong.\n')
      goThrough <- FALSE
    }
    
    return(goThrough)
  })






setMethod(
  f = "show",
  signature = "NG_graph",
  definition = function(object){
    cat("NG_graph object from ng_graph()\n")
    cat(paste("name:", object@name,"\n"))
    cat(paste("layout:", object@layout,"\n"))
  })


## layouting function
setMethod(
  f = "graphLayout",
  signature = "NG_graph",
  definition = function(graph,type,...){
    args <- list(...)
    
    ## random, kamadaKawaiSpring and fruchtermanReingold kamadaKawaiSpring and fruchtermanReingold currently not
    ## working for RBGL 1.32.1
    if(type %in% c("kamadaKawaiSpring", "fruchtermanReingold", "random")) {
      message(paste("Note: For graph \"", graph@name, "\": the ", type, " layout is temporarily not supported by the RBGL package, circle layout is used instead.", sep =""))
      type <- "circle"
    }
    
    ## ldist stands for label distance (from node)
    
    ## Distance from the Node Labels to the Node
    if(is(args$settings,"NG_Settings")){
      labelDistR <- args$settings@interaction@labelDistRadius
    }else if(is.null(args$ldist) == FALSE){
      labelDistR <- args$ldist
    }else{
      labelDistR <- 25
    }
    
    ## circle layout
    if(type == 'circle'){
      p <- numNodes(graph@graph)
      angle <- seq(0,2*pi, length=p+1)[1:p]
      
      ## radius of circle for nodes 
      r <- (min(graph@bbox) - 2*labelDistR - 2*graph@border)/2
      graph@xNodes <- cos(angle)*r + graph@bbox[1]/2
      graph@yNodes <- sin(angle)*r + graph@bbox[2]/2
      
      ## radius of circle for label names
      
      ##graph <- updateLabels(graph, ldist=labelDistR, orientation = 'center')
      
      r <- (min(graph@bbox) - 2*graph@border)/2
      graph@xLabels <- cos(angle)*r + graph@bbox[1]/2
      graph@yLabels <- sin(angle)*r + graph@bbox[2]/2
      
    }
    ## Depend on Rgraphviz
    ##			else if(any(type == c('neato','dot','twopi'))){
    ##							
    ##				## agopen is from the Rgraphviz library
    ##				vv <- agopen(graph = graph@graph, name = graph@name,
    ##						layout = type)
    ##				
    ##				x <- sapply(vv@AgNode,function(x)x@center@x)
    ##				y <- sapply(vv@AgNode,function(x)x@center@y)
    ##				nodeNames <- sapply(vv@AgNode,function(x)x@txtLabel@labelText)
    ##				
    ##				if(any(nodeNames != nodes(graph@graph))){
    ##					stop("[NG_graph:grahLayout] order of nodes after Layouting not the same anymore.\n")
    ##				}
    ##				
    ##				
    ##				## rescale graph
    ##				graph@xNodes <- (x-min(x))/(diff(range(x))) *
    ##						(graph@bbox[1]- 2*graph@border - 2*labelDistR)+graph@border+labelDistR
    ##				graph@yNodes <- (y-min(y))/(diff(range(y))) *
    ##						(graph@bbox[2]- 2*graph@border - 2*labelDistR)+graph@border+labelDistR
    ##				
    ##				graph <- updateLabels(graph, ldist=labelDistR, orientation = 'center')
    ##				
    ##				
    ##			}
    ## else if(type == 'random'){
    ##  
    ##   lay <- randomGraphLayout(graph@graph)
    ##  
    ##   if(all(colnames(lay) == nodes(graph@graph)) == FALSE){
    ##     stop("[NG_graph:graphLayout] random layout nodes not equal to layouted output.\n")
    ##   }
    ##  
    ##   x <- lay[1,]
    ##   y <- lay[2,]
    ##  
    ##   ##rescale graph
    ##   graph@xNodes <- (x-min(x))/(diff(range(x))) *
    ##     (graph@bbox[1]- 2*graph@border - 2*labelDistR)+graph@border+labelDistR
    ##   graph@yNodes <- (y-min(y))/(diff(range(y))) *
    ##     (graph@bbox[2]- 2*graph@border - 2*labelDistR)+graph@border+labelDistR
    ##  
    ##   ## graph@xNodes <- lay[1,]+graph@border
    ##   ## graph@yNodes <- lay[2,]+graph@border
    ##  
    ##   graph <- updateLabels(graph, ldist=labelDistR, orientation = 'center')
    ##  
    ## }
    ## Currently not working for RBGL 1.32.1
    ## else if(type == 'kamadaKawaiSpring'){
    ##  
    ##   lay <- kamadaKawaiSpringLayout( graph@graph, edge_or_side=1, es_length=1 )
    ##  
    ##   if(all(colnames(lay) == nodes(graph@graph)) == FALSE){
    ##     stop("[NG_graph:graphLayout] random layout nodes not equal to layouted output.\n")
    ##   }
    ##  
    ##   x <- lay[1,]
    ##   y <- lay[2,]
    ##   ##rescale graph
    ##   graph@xNodes <- (x-min(x))/(diff(range(x))) *
    ##     (graph@bbox[1]- 2*graph@border - 2*labelDistR)+graph@border+labelDistR
    ##   graph@yNodes <- (y-min(y))/(diff(range(y))) *
    ##     (graph@bbox[2]- 2*graph@border - 2*labelDistR)+graph@border+labelDistR
    ##  
    ##   graph <- updateLabels(graph, ldist=labelDistR, orientation = 'center')
    ##  
    ##  
    ## }else if(type == "fruchtermanReingold"){
    ##  
    ##   lay <- fruchtermanReingoldForceDirectedLayout(graph@graph)
    ##                                     # gursoyAtunLayout(graph@graph)$gursoyAtunLayout
    ##   if(all(colnames(lay) == nodes(graph@graph)) == FALSE){
    ##     stop("[NG_graph:graphLayout] random layout nodes not equal to layouted output.\n")
    ##   }
    ##  
    ##   x <- lay[1,]
    ##   y <- lay[2,]
    ##   ##rescale graph
    ##   graph@xNodes <- (x-min(x))/(diff(range(x))) *
    ##     (graph@bbox[1]- 2*graph@border - 2*labelDistR)+graph@border+labelDistR
    ##   graph@yNodes <- (y-min(y))/(diff(range(y))) *
    ##     (graph@bbox[2]- 2*graph@border - 2*labelDistR)+graph@border+labelDistR
    ##  
    ##   graph <- updateLabels(graph, ldist=labelDistR, orientation = 'center')
    ##  
    ##}
    else{
      stop(paste("[NG_graph:graphLayout] layout type \" ",type, "\" not known"))
    }
    
    
    
    graph@layout <- type
    return(graph)
  }
  )


## updateLabels
updateLabels <- function(graph, ldist, orientation = 'center'){
  if(is(graph,"NG_graph") == FALSE){
    stop("[NG_graph:updateLabels] graph is not from class NG_graph.\n")
  }
  
  if(orientation == 'same'){
    dx <- graph@xLabels - graph@xNodes 
    dy <- graph@yLabels - graph@yNodes
    
    dist <- apply(cbind(dx,dy),1,function(row)sqrt(sum(row^2)))
    
    ratio <- ldist/dist
    graph@xLabels <- graph@xNodes + ratio*dx
    graph@yLabels <- graph@yNodes + ratio*dy
    
    
  } else if(orientation == 'center'){
    center <- apply(cbind(graph@xNodes,graph@yNodes),2,mean)
    
    dx <- graph@xNodes-center[1]
    dy <- graph@yNodes-center[2]
    
    dist <- apply(cbind(dx,dy),1,function(row)sqrt(sum(row^2)))
    
    ratio <- (ldist+dist)/dist
    graph@xLabels <- center[1] + ratio*dx
    graph@yLabels <- center[2] + ratio*dy
  }
  
  return(graph)
}


## initialize
setMethod(
  f = "initialize",
  signature = "NG_graph",
  definition = function(.Object,name,graph,visualization,layout,sep){
    
    ## assign values
    .Object@name <- name
    if(is(graph,"graph")){
      .Object@graph <- graph
    }else{
      stop("[NG_graph:initialize] Graph is not from class \"graph\".")
    }
    .Object@bbox <- c(400,400)
    .Object@border <- 40
    
    
    .Object@layout <- layout
    
    .Object@sep <- sep
    
    
    ## Make graph Layout
    .Object <- graphLayout(.Object, layout)
    
    validObject(.Object)
    return(.Object)
  }
  )



## plot NG_graph
setMethod(
  f = "plot",
  signature = c(x = "NG_graph"),
  definition = function(x,y = NULL,...){
    plot(1,1, type = 'n', xlim = c(0,x@bbox[1]), ylim = c(0, x@bbox[2]),
         main = paste(x@name,' (\'',x@layout,'\' layout)', sep =''), axes = FALSE, xlab ='', ylab = '')
    box()
    axis(3,at = c(0,x@bbox[1]), labels = c(0,x@bbox[1]) )
    axis(2,at = c(0,x@bbox[2]), labels = c(x@bbox[2],0) )
    
    points(x@xNodes,x@yNodes*-1+x@bbox[2],pch = 19, cex = 2)
    text(x@xLabels,x@yLabels*-1+x@bbox[2],nodes(x@graph))
    
    
    edgeMatrix <- edgeMatrix(as(x@graph,"graphNEL"))
    
    segments(
      x@xNodes[edgeMatrix[1,]],
      x@yNodes[edgeMatrix[1,]]*-1+x@bbox[2],
      x@xNodes[edgeMatrix[2,]],
      x@yNodes[edgeMatrix[2,]]*-1+x@bbox[2])
    
  }
  )


## constructor function
ng_graph <- function(name, graph, sep=':', layout = 'circle'){
  
  if (grepl(' ', sep, fixed = TRUE)) {
    stop("[ng_graph] sep argument can not contain any spaces.")
  }
  
  if(!any(layout %in% c("circle", "kamadaKawaiSpring", "fruchtermanReingold", "random"))) {
    stop('[ng_graph] argument layout must be one of:\n \t "circle", "kamadaKawaiSpring", "fruchtermanReingold", "random"')
  }
  
  num <- sapply(gregexpr(sep, nodes(graph),fixed = TRUE),
                FUN = function(x){
                  v <- attr(x,'match.length')
                  if(all(length(v) ==1, v == -1)) {
                    0					
                  } else {
                    length(v)
                  }
                }
                )
  if (!all(num[1] == num)) {
    stop(paste("[ng_graph]: warning, the seperator number of occurences of sep in each graph node is\n",	
               paste(num, collapse = ", "),"\n"))
  }
  if(num[1] == 0){
    cat("[ng_graph]: warning, sep does not occur in some node names\n")
  }
  
  edgeDataDefaults(graph,"visited") <- FALSE
  return(new("NG_graph", name = name, graph = graph, layout = layout, sep = sep))
}





## resize the graph
setMethod(
  f = "resize",
  signature = c(object ="NG_graph", width = "numeric", height = "numeric"), 
  definition = function(object, width, height,...){
    args <- list(...)
    
    if(is.null(args$ldist)){
      ldist <- 25
    }else if(is(args$ldist,"NG_Settings")){
      ldist <- args$ldist@interaction@labelDistRadius
    }else{
      ldist <- args$ldist
    }
    
    ## adjust for border
    twoB <- 2*object@border
    
    object@xNodes <- (object@xNodes-object@border)/
      (object@bbox[1]-twoB)*(width-twoB)+object@border
    object@yNodes <- (object@yNodes-object@border)/
      (object@bbox[2]-twoB)*(height-twoB)+object@border
    
    object@xLabels <- (object@xLabels-object@border)/
      (object@bbox[1]-twoB)*(width-twoB)+object@border
    object@yLabels <- (object@yLabels-object@border)/
      (object@bbox[2]-twoB)*(height-twoB)+object@border
    
    object@bbox <- c(width,height)
    
    object <- updateLabels(object,ldist,orientation = 'same')
    
    return(object)
    
  })


setMethod(
  f = "adjacent",
  signature = "NG_graph",
  definition = function(graph,node,type = 'edge', retNr = TRUE){
    
    if(any(type == 'edge',type == 'node') == FALSE){
      stop('[ng_graph:adjacent] type argument wrong.\n')
    }
    
    if(length(node)!=1){
      stop('[ng_graph:adjacent] function needs only one node as an argument.\n')
    }
    
    nodeNames <- nodes(graph@graph)
    if(is.numeric(node)){
      nodeNr <- node
    }else{
      nodeNr <- unlist(match(node,nodeNames))
    }
    if(is.na(nodeNr)){
      stop('[ng_graph:adjacent] node not found.\n')
    }
    if(length(nodeNr)!=1){
      stop('[ng_graph:adjacent] node name not unique.\n')
    }
    
    
    adjNodesNames <- unlist(adj(graph@graph,node))
    if(all(type == 'node', retNr == FALSE)){
      return(adjNodesNames)
    }
    
    adjNodesNr <- match(adjNodesNames,nodeNames)
    if(all(type == 'node', retNr == TRUE)){
      return(adjNodesNr)
    }
    
    adjEdgesNr <- sapply(adjNodesNr,function(i){which(apply(graph@edgeMatrix,1,function(row){all(range(c(nodeNr,i)) == range(row))}))})
    if(all(type == 'edge', retNr == TRUE)){
      return(adjEdgesNr)
    }
    
    
    if(all(type == 'edge', retNr == FALSE)){
      return(edgeNames(graph@graph)[adjEdgesNr])
    }else{
      stop('[ng_graph:adjacent] check arguments.\n')
    }
  })




nodeNr <- function(graph,nodes){
  if(length(nodes)>1){
    nodesV <- nodes
  }else{
    nodesV <- unlist(strsplit(nodes, " "))
  }	
  match(nodesV, nodes(graph@graph))
  ## check whether each element is filled
}




setMethod(f = "ng_get",
          signature = "NG_graph",
          definition = function(obj, what=NULL, ...){
            possibleOptions <- c("name","graph","visitedEdges","layout")
            
            if(is.null(what)){
              cat("Get what? Possible options are: ")
              cat(paste(possibleOptions, collapse = ", "))
              cat("\n")
            }else{
              if(any(is.na(match(what,possibleOptions)))){
                stop(paste("[ng_get] object",what,"is not defined."))
              }
              
              if(what == "name"){
                return(obj@name)
              }else if(what == "graph"){
                return(obj@graph)
              }else if(what == "data"){
                return(obj@data)
              }else if(what == "layout"){
                return(obj@layout)
              }
            }
          })


setMethod(f = "ng_set",
          signature = "NG_graph",
          definition = function(object){
            possibleOptions <- c("name","graph","layout")
            paste("Replace what? Possible options are:", paste(possibleOptions, collapse = ", "),'\n')
            cat('Use ng_set<- to set a value.\n')
          })


setReplaceMethod(
  f = "ng_set","NG_graph",
  function(object,what,value){
    possibleOptions <- c("name","graph","layout")
    tmp <- object
    
    if(!(what %in% possibleOptions)) {
      stop(paste("Replace what? Possible options are: ", paste(possibleOptions, collapse = ", "),'\n'))				
    }	
    
    slot(object, what, check = TRUE) <- value 
    if(validObject(object)){
      return(object)
    }else {
      cat("[ng_set] assignment is wrong\n")
      return(tmp)
    } 
  })
