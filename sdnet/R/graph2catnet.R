#########################################################################
# import network structure from a graph object

edges2catnet <- function(.Object, nodes, edges, maxcats=2, cats=NULL, probs=NULL) {

  ## Nodes and Parents  
  .Object@nodes <- nodes
  .Object@meta <- ""
  nn <- length(.Object@nodes)
  .Object@numnodes <- nn
  pars <- vector("list", nn)  

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
        if(length(pars) < k)  
          pars[[k]] <- c(i)  
        else  
          pars[[k]] <- c(pars[[k]], i) 
      } 
      else {  
        stop("Invalid node in edge (", i, ", ", j, ")\n")  
      }  
    }  
  }

  maxpars <- 0  
  for(i in (1:length(pars))) {  
    if(maxpars < length(pars[[i]]))  
      maxpars <- length(pars[[i]])  
  }  
  .Object@maxpars <- as.integer(maxpars)    
  .Object@pars <- pars  

  ## Categories  
  .Object@cats <- list(NULL)  
  if(length(cats) == nn) {  
    maxcats <- 0  
    for(i in (1:nn)) {  
      if(maxcats < length(cats[i]))  
        maxcats <- length(cats[i])  
      .Object.cats[[i]] <- cats[i]  
    }  
  }  
  else{
    if(maxcats < 2) {
      maxcats <- 2
      warning("Set maxcats to 2")
    }
    for(i in (1:nn)) {  
      .Object@cats[[i]] <- as.character(1:maxcats)
    }  
  }  
  .Object@maxcats <- as.integer(maxcats)  
  ## Probability Matrix  
  .Object@probs <- list(NULL)  
  if(!is.null(probs)) {  
    .Object@probs <- probs  
  }  
  else {              
    #poutlist <- NULL
    #poutlist <- lapply(seq(1,nn), function(parid) {
    #  setRandomProb(parid, .Object@pars[[parid]], .Object@cats,
    #                seq(1, length(.Object@pars[[parid]])))
    #})
    #.Object@probs <- poutlist
  }

  .Object@nodecomplx <- sapply(1:.Object@numnodes, function(x) nodeComplexity(.Object, x))
  .Object@complx <- as.integer(sum(.Object@nodecomplx))
  
  .Object@loglik <- 0
  .Object@nodelik <- NA
  
  return(.Object)  
}

# inherit catnet from graph
graph2catnet <- function(object, gr, maxcats=2, cats=NULL, probs=NULL) {  
  if(!isS4(gr))
    stop("Not a valid graph")
  id <- which(slotNames(gr) == "nodes")
  if(length(id) < 1)
    stop("No slot names in graph")
  id <- which(slotNames(gr) == "edgeL")
  if(length(id) < 1)
    stop("No slot edgeL in graph")
  return(edges2catnet(object, gr@nodes, gr@edgeL, maxcats, cats, probs))
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


## generate a DAG graph
## returns a graphNEL object
#genRandomGraph <- function(numnodes, maxparents) {
#  idx <- sample(seq(1, numnodes))
#  pars <- vector("list", numnodes)
#  for(i in (2:numnodes)) {
#    npars <- floor(maxparents*runif(1,0,1) + 0.5)
#    if(npars == 0)
#      next;
#    if(npars > i-1)
#      npars <- i - 1
#    pars[[idx[i]]] <- sample(idx[1:(i-1)], npars)
#  }
#  nodes <- sapply(seq(1,numnodes), function(i) paste("N", i ,sep=""))
#  edges <- vector("list", numnodes)
#  for(i in (1:numnodes)) {
#    for(j in (1:numnodes)) {
#       if(is.null(pars[[j]]))
#         next
#       ll <- which(pars[[j]]==i)
#       if(is.na(ll[1]))
#         next
#       edges[[i]] <- c(edges[[i]], nodes[j])
#     }
#  }
#  edges<-setNames(edges, rep("edges", numnodes))
#  edgesL <- sapply(1:numnodes, function(j, edges) {
#    list(edges[j])
#  }, edges)
#  edgesL <- setNames(edgesL, nodes)
#  if(require("graph"))
#    return(try(new("graphNEL", nodes=nodes, edgeL=edgesL, edgemode="directed")))
#  else
#    return(NULL)
#}

listGraphEdges <- function(object) {
  if(!is(object, "catNetwork") || object@numnodes < 1)
    return(NULL)
  numnodes <- object@numnodes
  nodes <- object@nodes
  pars <- object@pars
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

