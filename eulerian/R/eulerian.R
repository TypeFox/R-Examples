require(graph)

eulerian <- function(graph, start=NULL){
  
  if(!inherits(graph,"graphNEL"))
    stop("'graph' must be a graphNEL object.")
  
  nodes <- nodes(graph)
  
  if(!is.null(start)){
    if(!inherits(start, "character") || length(start)!=1)
      stop("start must be NULL or a character of length 1.")
    if(!start %in% nodes)
      stop("start must be either NULL or a node of graph.")
  }
  
  if(!hasEulerianPath(graph, start)){
    errmsg <- ifelse(is.null(start), "There is no eulerian path.", paste0("There is no eulerian path starting from ",start,"."))
    stop(errmsg)
  }
  
  # if no nodes/edges in graph, return an empty path
  if(numNodes(graph)==0 || (numEdges(graph)==0 && is.null(start)))
    return(character())
  
  # algorithm depends on directed/undirected graph
  directed <- isDirected(graph)
  
  # remove nodes without any edge
  deg <- degree(graph, Nodes=nodes)
  if(directed){
    nodesToRemove <- nodes[(deg$outDegree + deg$inDegree)==0]
  } else{
    nodesToRemove <- nodes[deg==0]
  }
  
  if(length(nodesToRemove)>0){
    graph <- removeNode(object=graph, node=nodesToRemove)
    nodes <- nodes(graph)
  }
  
  
  selectAStartNode <- function(graph){
    startNode <- NA
    nodes <- nodes(graph)
    deg <- degree(graph, Nodes=nodes)
    
    if(isDirected(graph)){
      degDiff <- deg$outDegree - deg$inDegree
      degDiff1Index <- which(degDiff==1)
      if(length(degDiff1Index)==1){
        startNode <- nodes[degDiff1Index]
      }else if(length(degDiff1Index)==0){
        startNode <- nodes[1]
      }  
    } else{
      oddDegIndex <- which(deg%%2==1)
      if(length(oddDegIndex)==0){
        startNode <- nodes[1]
      } else if(length(oddDegIndex) <= 2){
        startNode <- nodes[oddDegIndex[1]]
      }
    }
    
    return(startNode)
  }
  
  eulerPath <- c()
  ed <- edges(graph)
  startIndex <- 1
  endIndex <- length(ed)
  curNode <- ifelse(is.null(start), selectAStartNode(graph), start)
  
  while(!is.na(curNode)){
    epath <- curNode
    
    while(TRUE){
      # Random walk
      moves <- ed[[curNode]]
      nextNode <- ifelse(length(moves)>0, moves[1], NA)
      if(is.na(nextNode))
        break;
      
      # remove edge from curNode to nextNode
      ed[[curNode]] <- ed[[curNode]][-1]
      if(!directed && curNode!=nextNode){
        removeIndex <- which(ed[[nextNode]]==curNode)[1]
        ed[[nextNode]] <- ed[[nextNode]][-removeIndex]
      }
      epath <- append(epath, nextNode)
      curNode <- nextNode
    }
    
    # update euler path
    if(length(eulerPath)==0){
      eulerPath <- epath
    }
    else{
      insertIndex <- which(eulerPath==epath[1])[1]
      eulerPath <- append(eulerPath,after=insertIndex,values=epath[-1])
    }
    
    #select next start node (still unexplored)
    curNode <- NA
    while((startIndex <- startIndex + 1) <= endIndex){
      node <- eulerPath[startIndex]
      if(length(ed[[node]])>0){
        curNode <- node
        break
      }
    }
    
  }
  
  return(eulerPath)
}

hasEulerianPath <- function(graph, start=NULL){
  
  if(!inherits(graph,"graphNEL"))
    stop("'graph' must be a graphNEL object.")
  
  if(!is.null(start)){
    if(!inherits(start, "character") || length(start)!=1)
      stop("start must be NULL or a character of length 1.")
    if(!start %in% nodes(graph))
      stop("start must be either NULL or a node of graph.")
  }
  
  nodes <- nodes(graph)
  if(length(nodes)==0)
    return(TRUE)
  if(numEdges(graph)==0 && is.null(start))
    return(TRUE)
  
  # algorithm varies for directed and undirected graphs
  directed <- isDirected(graph)
  
  
  # the nodes, that have at least 1 edge, must be connected
  deg <- degree(graph, Nodes=nodes)
  if(directed){
    nodes <- nodes[(deg$outDegree + deg$inDegree)>0]
  } else{
    nodes <- nodes[deg>0]
  }
  graph <- subGraph(graph=graph, snodes=nodes)
  if(!isConnected(graph))
    return(FALSE)
  
  hasEulerPath <- FALSE
  if(directed){
    # each nodes' indegree and outdegree must be same (with possible exception at start & end node)
    deg <- degree(graph, Nodes=nodes)
    degDiff <- deg$outDegree - deg$inDegree
    degDiff1Index <- which(degDiff==1)
    if(length(degDiff1Index)==1){
      startNode <- nodes[degDiff1Index]
      if(is.null(start) || start==startNode){
        hasEulerPath <- TRUE
      }
    }else if(length(degDiff1Index)==0){
      if(is.null(start) || start %in% nodes){
        hasEulerPath <- TRUE
      }
    }
  } else {
    # each nodes' degree must be even (with possible exception at start & end node)
    deg <- degree(graph, Nodes=nodes)
    edL <- edges(graph, which=nodes)
    selfDeg <- sapply(nodes, function(n){
      return(length(which(edL[[n]]==n)))
    })
    deg <- deg + selfDeg
    oddDegIndex <- which(deg%%2==1)
    if(length(oddDegIndex)==0){
      if(is.null(start) || start %in% nodes){
        hasEulerPath <- TRUE
      }
    } else if(length(oddDegIndex) <= 2){
      startNodes <- nodes[oddDegIndex]
      if(is.null(start) || start %in% startNodes){
        hasEulerPath <- TRUE
      }
    }
  }
  
  return(hasEulerPath)
}

hasEulerianCycle <- function(graph){
  
  if(!inherits(graph,"graphNEL"))
    stop("'graph' must be a graphNEL object.")
  
  nodes <- nodes(graph)
  if(length(nodes)==0)
    return(TRUE)
  if(numEdges(graph)==0)
    return(TRUE)
  
  # algorithm varies for directed and undirected graphs
  directed <- isDirected(graph)
  
  # the nodes, that have at least 1 edge, must be connected
  deg <- degree(graph, Nodes=nodes)
  if(directed){
    nodes <- nodes[(deg$outDegree + deg$inDegree)>0]
  } else{
    nodes <- nodes[deg>0]
  }
  graph <- subGraph(graph=graph, snodes=nodes)
  if(!isConnected(graph))
    return(FALSE)
  
  hasEulerCycle <- FALSE
  if(directed){
    # every nodes' indegree and outdegree must be same
    deg <- degree(graph, Nodes=nodes)
    degDiff <- deg$outDegree - deg$inDegree
    hasEulerCycle <- all(degDiff==0)
  } else {
    # every nodes' degree must be even
    deg <- degree(graph, Nodes=nodes)
    edL <- edges(graph, which=nodes)
    selfDeg <- sapply(nodes, function(n){
      return(length(which(edL[[n]]==n)))
    })
    deg <- deg + selfDeg
    hasEulerCycle <- all(deg%%2==0)
  }
  
  return(hasEulerCycle)
}
