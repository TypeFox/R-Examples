#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2013 the statnet development team
######################################################################
# checks the network attributes to see if a vertex.pid nodal attribute is specified
vertex.pid.check <- function(nd) {
  vpid.name <- get.network.attribute(nd, 'vertex.pid')
  if (is.null(vpid.name)){
    warning("Input network does not have a 'vertex.pid' attribute giving the name of the vertex attribute to be used as a vertex.pid")
    return(FALSE)
  }
  
  if (vpid.name %in% list.vertex.attributes(nd)) {
    vertex.pid.vals <- get.vertex.attribute(nd, vpid.name)
    if (length(unique(vertex.pid.vals)) != network.size(nd)){
      stop("vertex.pid attribute must be specified and unique for each vertex.")
    }
    if (any(is.na(vertex.pid.vals))){
      stop("vertex.pid attribute must not be NA for any vertex.")
    }
    
  } else {
    stop(paste("Missing vertex.pids: Each vertex must have a vertex.pid attribute named '",vpid.name,"'"))
  }
  
  return(TRUE)
}

# checks the network attributes to see if a edge.pid nodal attribute is specified
edge.pid.check <- function(nd) {
  epid.name <- get.network.attribute(nd, 'edge.pid')
  if (is.null(epid.name)){
    warning("Input network does not have an 'edge.pid' attribute giving the name of the edge attribute to be used as an edge.pid")
    return(FALSE)
  }
  
  if (epid.name %in% list.edge.attributes(nd)) {
    edge.pid.vals <- get.edge.value(nd, epid.name)
    if (length(unique(edge.pid.vals)) != network.edgecount(nd)){
      stop("edge.pid attribute must be specified and unique for each edge.")
    }
    if (any(is.null(edge.pid.vals))){
      stop("edge.pid attribute must not be NULL for any existing edge")
    }
    if (any(is.na(edge.pid.vals))){
      stop("edge.pid attribute must not be NA for any edge")
    }
    
  } else {
    stop(paste("Missing edge.pids: Each edge must have a edge.pid attribute named '",epid.name,"'"))
  }
  return(TRUE)
}

# takes a networkDynamic object and a vertex.pid, returns the vertex id (number) of the node
# returns NA if not found
get.vertex.id <- function(nd, pid) {
  pid.name <- get.network.attribute(nd, 'vertex.pid')
  if (is.null(pid.name)){
    stop("Input network does not have a 'vertex.pid' attribute giving the name of the vertex attribute to be used as a vertex.pid")
  # pid.name = 'vertex.pid'
  }
  
  res <- match(pid,get.vertex.attribute(nd, pid.name))
  if (length(res)==0) return(NA)
  res
}


# takes a networkDynamic object and a vertex id, returns the vertex pid attribute for the node
# returns NA if not found
get.vertex.pid <- function(nd, id) {
  pid.name <- get.network.attribute(nd, 'vertex.pid')
  if (is.null(pid.name)){
    # pid.name <- 'vertex.pid'
    stop("Input network does not have a 'vertex.pid' attribute giving the name of the vertex attribute to be used as a vertex.pid")
  } 
  get.vertex.attribute(nd, pid.name)[id]
}

# takes a networkDynamic object and an edge.pid, returns the edgex id (number) of the edge
# returns NA if not found
get.edge.id <- function(nd, pid) {
  pid.name <- get.network.attribute(nd, 'edge.pid')
  if (is.null(pid.name)){
    stop("Input network does not have an 'edge.pid' attribute giving the name of the edge attribute to be used as a edge.pid")
  }
  ndPids<-get.edge.value(nd, pid.name,unlist=FALSE)
  res <- sapply(pid,match,ndPids,USE.NAMES=FALSE)
  if (length(res)==0) return(NA)
  res
}

# takes a networkDynamic object and a edge id, returns the edge pid attribute for the node
# returns NA if not found
get.edge.pid <- function(nd, id) {
  pid.name <- get.network.attribute(nd, 'edge.pid')
  if (is.null(pid.name)){
    stop("Input network does not have an 'edge.pid' attribute giving the name of the edge attribute to be used as an edge.pid")
  } 
  sapply(get.edge.value(nd, pid.name,unlist=FALSE)[id],function(x){ifelse(is.null(x),NA,x)})
}

# adding vertex.pid to new node if specified
# this function needs to overide the version in network to make sure
# that the pid is hanled correctly
add.vertices.networkDynamic <- function(x, nv, vattr=NULL, last.mode=TRUE, vertex.pid=NULL, ...) {
  
  # todo: test x is net
  
  # for modify in place
  xn <- substitute(x)
  
  if (last.mode || !is.bipartite(x)) {
    n <- get.network.attribute(x, "n")
  } else {
    # adding to the first mode in bipartite
    # "Where bipartite > 0, network methods will
    #  automatically assume that vertices with IDs less than or equal to bipartite belong
    #  to one such class, with those with IDs greater than bipartite belonging to the other."
    n <- x %n% "bipartite"
  }
  newids <- (n+1):(n+nv)
  
  # only mess with vertex pids if it is named at the network level
  pid.name <- get.network.attribute(x, 'vertex.pid')
  if (!is.null(pid.name)){
    if( !is.null(vertex.pid)) {
      pid.name <- get.network.attribute(x, 'vertex.pid')
      if(length(unique(vertex.pid))!=length(vertex.pid)){
        stop("All vertex.pid values must be unique")
      }
      
      if (length(vertex.pid) != nv){
        stop("Number of Vertex.pids supplied does not match number of new vertices added")
      }
    } else {  #using pids, but not given so need to make some up
      vertex.pid<-generatePids(nv,get.vertex.attribute(x,pid.name))
    }
  }
  
  x<-add.vertices.network(x=x, nv=nv, vattr=vattr, last.mode=last.mode)
  # TODO: if there is an error condition below, network will allready have been modified. 
  
  # only mess with vertex pids if it is named at the network level
  if (!is.null(pid.name)){
    set.vertex.attribute(x, pid.name, vertex.pid, v=newids)
    vertex.pid.check(x)
  }
  
  # modify-in-place voodo
  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, x)))
  invisible(x)
}

# override add.edges in network to manage pid stuff
add.edges.networkDynamic<-function(x, tail, head, names.eval=NULL, vals.eval=NULL,edge.pid=NULL,...){
  # for modify in place
  xn <- substitute(x)
  
  mnextOrig<-x%n%'mnext'
  x<-add.edges.network(x=x,tail=tail,head=head,names.eval=names.eval,vals.eval=vals.eval,...=...)
  mnextNew<-x%n%'mnext'
  
  # only mess with edge pids if it is named at the network level
  pid.name <- get.network.attribute(x, 'edge.pid')
  if (!is.null(pid.name)){
    if( !is.null(edge.pid)) {
      pid.name <- get.network.attribute(x, 'edge.pid')
      if(length(unique(edge.pid))!=length(edge.pid)){
        stop("All edge.pid values must be unique")
      }
      
      if (length(edge.pid) != length(tail)){
        stop("Number of edge.pids supplied does not match number of new edges added")
      }
    } else {  #using pids, but not given so need to make some up
      edge.pid<-generatePids(length(tail),get.edge.value(x,pid.name))
    }
  }
  
  # only mess with vertex pids if it is named at the network level
  if (!is.null(pid.name)){
    set.edge.attribute(x, pid.name, edge.pid, e=seq.int(from=mnextOrig,length.out=(mnextNew-mnextOrig)))
    edge.pid.check(x)
  }
  
  # modify-in-place voodo
  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, x)))
  invisible(x)
}

# override add.edges in network to manage pid stuff
add.edge.networkDynamic<-function(x, tail, head, names.eval=NULL, vals.eval=NULL,edge.check=FALSE,edge.pid=NULL,...){
  # for modify in place
  xn <- substitute(x)
  ev <- parent.frame()
  
  mnextOrig<-x%n%'mnext'
  x<-add.edge.network(x=x,tail=tail,head=head,names.eval=names.eval,vals.eval=vals.eval,edge.check=edge.check,...=...)
  mnextNew<-x%n%'mnext'
  
  # only mess with edge pids if it is named at the network level
  pid.name <- get.network.attribute(x, 'edge.pid')
  if (!is.null(pid.name)){
    if( !is.null(edge.pid)) {
      pid.name <- get.network.attribute(x, 'edge.pid')
      if (length(edge.pid) >1){
        stop("Only one edge.pid can be specified for add.edge")
      }
    } else {  #using pids, but not given so need to make some up
      edge.pid<-generatePids(1,get.edge.value(x,pid.name))
    }
  }
  
  # only mess with vertex pids if it is named at the network level
  if (!is.null(pid.name)){
    set.edge.attribute(x, pid.name, edge.pid, e=seq.int(from=mnextOrig,length.out=(mnextNew-mnextOrig)))
    edge.pid.check(x)
  }
  
  # modify-in-place voodo
  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, x)))
  invisible(x)
}

# # don't yet have a magic way to do this, trying generating a random hex string and testing that it is unique. Hopefully randomness assumes they will never be reused, and the tempfile function is robust enough.  could always overide this function if want to generate in a new way. 
generatePids<- function(n,existingIds=character(0)){
  if(anyDuplicated(existingIds)){
    stop(".generatePids cannot generate a new unique pid because existing pids are not unique")
  }
  newids <-replicate(n,substring(tempfile(pattern='',tmpdir=''),first=2))
  while(anyDuplicated(c(newids,existingIds))){
    newids <-replicate(n,substring(tempfile(pattern='',tmpdir=''),first=2))
  }
  return(newids)
}

# create sets of pids for vertices and edges
initialize.pids<-function(nd){
  if(!is.network(nd)){
    stop("Persistant ids can not be added to objects that are not networks")
  }
  # for modify in place
  xn <- substitute(nd)
  
  if (is.null(nd%n%'vertex.pid')){
    vpids<-generatePids(network.size(nd))
    set.network.attribute(nd,'vertex.pid','vertex.pid')
    set.vertex.attribute(nd,'vertex.pid',vpids)
  } else {
    warning("Network already contains a vertex.pid attribute, no new persistant ids were generated for vertices")
  }
  if (is.null(nd%n%'edge.pid')){
    epids<-generatePids(network.edgecount(nd))
    set.network.attribute(nd,'edge.pid','edge.pid')
    set.edge.attribute(nd,'edge.pid',epids)
  } else {
    warning("Network already contains a edge.pid attribute, no new persistant ids were generated for edges")
  }
  # modify-in-place voodo
  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, nd)))
  invisible(nd)
}
