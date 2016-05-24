#########################################################################
# import network structure from a graph object

edges2catnet <- function(.Object, nodes, edges, maxcats=2, cats=NULL, probs=NULL) {

  ## Nodes and Parents  
  .Object@nodes <- nodes
  .Object@meta <- ""
  nn <- length(.Object@nodes)
  .Object@numnodes <- nn
  parents <- vector("list", nn)  

  for(i in (1:nn)) {
    nedges <- edges[[i]]
    if(length(nedges) == 0)   
      next
    if(is.list(nedges))
      nedges <- nedges[[1]]
    if(length(nedges) == 0)   
      next
    for(j in (1:length(nedges))) {  
      k <- which(.Object@nodes == nedges[j])  
      if(length(k)==1) {  
        if(length(parents) < k)  
          parents[[k]] <- c(i)  
        else  
          parents[[k]] <- c(parents[[k]], i) 
      } 
      else {  
        stop("Invalid node in edge (", i, ", ", j, ")\n")  
      }  
    }  
  }

  maxpars <- 0  
  for(i in (1:length(parents))) {  
    if(maxpars < length(parents[[i]]))  
      maxpars <- length(parents[[i]])  
  }  
  .Object@maxParents <- as.integer(maxpars)    
  .Object@parents <- parents  

  ## Categories  
  .Object@categories <- list(NULL)  
  if(length(cats) == nn) {  
    maxcats <- 0  
    for(i in (1:nn)) {  
      if(maxcats < length(cats[i]))  
        maxcats <- length(cats[i])  
      .Object.categories[[i]] <- cats[i]  
    }  
  }  
  else{
    if(maxcats < 2) {
      maxcats <- 2
      warning("Set maxCategories to 2")
    }
    for(i in (1:nn)) {  
      .Object@categories[[i]] <- as.character(1:maxcats)
    }  
  }  
  .Object@maxCategories <- as.integer(maxcats)  
  ## Probability Matrix  
  pp <- list(NULL)  
  if(!is.null(probs)) {  
    pp <- probs  
  }  
  else {              
    poutlist <- NULL
    poutlist <- lapply(seq(1,nn), function(parid) {
      setRandomProb(parid, .Object@parents[[parid]], .Object@categories,
                    seq(1, length(.Object@parents[[parid]])))
    })
    .Object@probabilities <- poutlist
  }

  .Object@nodeComplexity <- sapply(1:.Object@numnodes, function(x) nodeComplexity(.Object, x))
  .Object@complexity <- as.integer(sum(.Object@nodeComplexity))
  
  .Object@likelihood <- 0
  .Object@nodeLikelihood <- NA
  
  return(.Object)  
}

# inherit catnet from graph
graph2catnet <- function(object, gr, maxCategories=2, cats=NULL, probs=NULL) {  
  if(!isS4(gr))
    stop("Not a valid graph")
  id <- which(slotNames(gr) == "nodes")
  if(length(id) < 1)
    stop("No slot names in graph")
  id <- which(slotNames(gr) == "edgeL")
  if(length(id) < 1)
    stop("No slot edgeL in graph")
  return(edges2catnet(object, gr@nodes, gr@edgeL, maxCategories, cats, probs))
}

setRandomProbMatrixForm <- function(idroot, ppars, pcatlist, idx, poutlist) {
  if(is.null(ppars) || length(idx) < 1) {
    if(length(pcatlist[[idroot]]) < 2) {
      poutlist <- c(poutlist, 1)
      return(NULL)
    }
    plist <- sapply(pcatlist[[idroot]], function(x) runif(1,0.01,0.99))
    plist <- plist/sum(plist)
    plist <- sapply(seq(1,length(plist)), function(n, plist){
      plist[n] <- floor(100*plist[n])/100
      }, plist)
    plist[1] <- 1 - sum(plist[-1])
    return(plist)
  }
  else {
    id <- ppars[idx[1]]
    poutlist <- sapply(pcatlist[[id]],
           function(cat, idroot1, ppars1, pcatlis1t, idx1, poutlist1)
           setRandomProbMatrixForm(idroot1, ppars1, pcatlis1t, idx1, poutlist1),
           idroot, ppars, pcatlist, idx[-1], poutlist)
  }
}


listGraphEdges <- function(object) {
  if(!is(object, "catNetwork") || object@numnodes < 1)
    return(NULL)
  numnodes <- object@numnodes
  nodes <- object@nodes
  pars <- object@parents
  edges <- vector("list", numnodes)
  for(i in (1:numnodes)) {
    for(j in (1:numnodes)) {
      if(is.null(pars[[j]]))
        next
      idx <- which(pars[[j]]==i)
      if(is.na(idx[1]))
        next
      edges[[i]] <- c(edges[[i]], j)
    }
  }
  edges<-setNames(edges, rep("edges", numnodes))
  edgesL <- sapply(1:numnodes, function(j, edges) {
    list(edges[j])
  }, edges)
  edgesL <- setNames(edgesL, nodes)
  return(edgesL)
}

