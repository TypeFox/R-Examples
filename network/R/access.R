######################################################################
#
# access.R
#
# Written by Carter T. Butts <buttsc@uci.edu>; portions contributed by
# David Hunter <dhunter@stat.psu.edu> and Mark S. Handcock
# <handcock@u.washington.edu>.
#
# Last Modified 02/26/13
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/network package
#
# This file contains various routines for accessing network class objects.
#
# Contents:
#
#   add.edge
#   add.edges
#   add.vertices
#   delete.edge.attribute
#   delete.edges
#   delete.network.attribute
#   delete.vertex.attribute
#   delete.vertices
#   get.edge.attribute
#   get.edge.value
#   get.edgeIDs
#   get.edges
#   get.inducedSubgraph
#   get.network.attribute
#   get.neighborhood
#   get.vertex.attribute
#   has.loops
#   is.adjacent
#   is.bipartite
#   is.directed
#   is.hyper
#   is.multiplex
#   is.network
#   list.edge.attributes
#   list.network.attributes
#   list.vertex.attributes
#   network.dyadcount
#   network.edgecount
#   network.naedgecount
#   network.size
#   network.vertex.names
#   network.vertex.names<-
#   permute.vertexIDs
#   set.edge.attribute
#   set.edge.value
#   set.network.attribute
#   set.vertex.attribute
#   valid.eids
#
######################################################################


#Add a single edge to a network object.
# S3 method dispatch for add edge
add.edge<-function(x, tail, head, names.eval=NULL, vals.eval=NULL, edge.check=FALSE, ...){
  xn<-substitute(x)
  UseMethod("add.edge") 
  if(.validLHS(xn,parent.frame())){  #If x not anonymous, set in calling env 
    on.exit(eval.parent(call('<-',xn,x)))
  }
  invisible(x) 
} 

add.edge.network<-function(x, tail, head, names.eval=NULL, vals.eval=NULL, edge.check=FALSE, ...){ 
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("add.edge requires an argument of class network.")
  #Do the deed
  xn<-substitute(x)
  if(.validLHS(xn,parent.frame())){  #If x not anonymous, set in calling env 
    on.exit(eval.parent(call('<-',xn,x)))
  }
  x<-.Call(addEdge_R,x,tail,head,names.eval,vals.eval,edge.check)
  invisible(x)
}

# S3 method dispatch for add.edges
add.edges<-function(x, tail, head, names.eval=NULL, vals.eval=NULL, ...){
  xn<-substitute(x)
  UseMethod("add.edges") 
  if(.validLHS(xn,parent.frame())){  #If x not anonymous, set in calling env 
    on.exit(eval.parent(call('<-',xn,x)))
  }
  invisible(x) 
} 


# Add multiple edges to network x.  Tail must be a list, each element of
# which is the tail set for a given edge (ditto for head).  If edge values
# are provided, they must be given similarly as lists of lists.
add.edges.network<-function(x, tail, head, names.eval=NULL, vals.eval=NULL, ...){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("add.edges requires an argument of class network.")
  #Ensure that the inputs are set up appropriately 
  if(!is.list(tail))
    tail<-as.list(tail)
  if(!is.list(head))
    head<-as.list(rep(head,length=length(tail)))
  if(is.null(names.eval))
    names.eval<-replicate(length(tail),NULL)
  else if(!is.list(names.eval))
    names.eval<-as.list(rep(names.eval,length=length(tail)))
  if(is.null(vals.eval))
    vals.eval<-replicate(length(tail),NULL)
  else if(!is.list(vals.eval))
    vals.eval<-as.list(rep(vals.eval,length=length(names.eval)))
  if(length(unique(c(length(tail),length(head),length(names.eval), length(vals.eval))))>1)
    stop("head, tail, names.eval and vals.eval lists passed to add.edges must be of the same length!\n")
  edge.check<-list(...)$edge.check
  if(is.null(edge.check))
    edge.check<-FALSE
  #Pass the inputs to the C side
  xn<-substitute(x)
  x<-.Call(addEdges_R,x,tail,head,names.eval,vals.eval,edge.check)
  if(.validLHS(xn,parent.frame())){  #If x not anonymous, set in calling env 
    on.exit(eval.parent(call('<-',xn,x)))
  }
  invisible(x)
}



# S3 method dispatch for add.vertices
add.vertices<-function(x, nv, vattr=NULL, last.mode=TRUE, ...){
  xn<-substitute(x)
  UseMethod("add.vertices") 
  if(.validLHS(xn,parent.frame())){  #If x not anonymous, set in calling env 
    on.exit(eval.parent(call('<-',xn,x)))
  }
  invisible(x) 
}

# Add nv vertices to network x.  Vertex attributes (in addition to those which
# are required) are to be provided in vattr; vattr must be a list containing
# nv elements, each of which is equal to the desired val[i] entry.
add.vertices.network<-function(x, nv, vattr=NULL, last.mode=TRUE, ...){ 
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("add.vertices requires an argument of class network.\n")
  #Check the vertex attributes, to be sure that they are legal
  if(!is.null(vattr)){
    if(is.list(vattr))
      vattr<-rep(vattr,length=nv)
    else
      vattr<-as.list(rep(vattr,length=nv))
  }
  #Perform the addition
  xn<-substitute(x)
  if(nv>0){
    if(.validLHS(xn,parent.frame())){  #If x not anonymous, set in calling env 
      on.exit(eval.parent(call('<-',xn,x)))
    }
    if(last.mode||(!is.bipartite(x))){
      x<-.Call(addVertices_R,x,nv,vattr)
    }else{
      
      nr<-nv
      nc<-0
      nnew<-nr+nc
      nold<-network.size(x)
      bip<-x%n%"bipartite"
      x<-.Call(addVertices_R, x, nv, vattr)
      
      if(nr>0){
        if(bip>0)
          orow<-1:bip
        else
          orow<-NULL
        if(nold-bip>0)
          ocol<-(bip+1):nold
        else
          ocol<-NULL
        
        ncol<-NULL
        nrow<-(nold+nnew-nr+1):(nold+nnew)
        permute.vertexIDs(x,c(orow,nrow,ocol,ncol))
        set.network.attribute(x,"bipartite",bip+nr)
      }
    }
  }

  invisible(x)
}


# Remove all instances of the specified attribute(s) from the edge set
#
delete.edge.attribute<-function(x,attrname){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("delete.edge.attribute requires an argument of class network.")
  #Remove the edges
  xn<-substitute(x)
  x<-.Call(deleteEdgeAttribute_R,x,attrname)
  if(.validLHS(xn,parent.frame())){  #If x not anonymous, set in calling env 
    on.exit(eval.parent(call('<-',xn,x)))
  }
  invisible(x)
}
 
 
# Remove specified edges from the network.
#
delete.edges<-function(x,eid){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("delete.edges requires an argument of class network.")
  xn<-substitute(x)
  if(length(eid)>0){
    #Perform a sanity check
    if((min(eid)<1)|(max(eid)>length(x$mel)))
      stop("Illegal edge in delete.edges.\n")
    #Remove the edges
    x<-.Call(deleteEdges_R,x,eid)
    if(.validLHS(xn,parent.frame())){  #If x not anonymous, set in calling env 
      on.exit(eval.parent(call('<-',xn,x)))
    }
  }
  invisible(x)
}

# Remove the specified network-level attribute(s)
#
delete.network.attribute<-function(x,attrname){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("delete.network.attribute requires an argument of class network.")
  #Remove the edges
  xn<-substitute(x)
  x<-.Call(deleteNetworkAttribute_R,x,attrname)
  if(.validLHS(xn,parent.frame())){  #If x not anonymous, set in calling env 
    on.exit(eval.parent(call('<-',xn,x)))
  }
  invisible(x)
}


# Remove all instances of the specified attribute(s) from the vertex set
#
delete.vertex.attribute<-function(x,attrname){
  #Check to be sure we were called with a network, and that it has vertices
  if(!is.network(x))
    stop("delete.vertex.attribute requires an argument of class network.")
  #Remove the attribute (or do nothing, if there are no vertices)
  if(network.size(x)>0){
    xn<-substitute(x)
    x<-.Call(deleteVertexAttribute_R,x,attrname)
    if(.validLHS(xn,parent.frame())){  #If x not anonymous, set in calling env 
      on.exit(eval.parent(call('<-',xn,x)))
    }
  }
  invisible(x)
}


# Remove specified vertices (and associated edges) from the network.
#
delete.vertices<-function(x,vid){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("delete.vertices requires an argument of class network.")
  #Remove any vids which are out of bounds
  vid<-vid[(vid>0)&(vid<=network.size(x))]
  #Do the deed, if still needed
  xn<-substitute(x)
  if(length(vid)>0){
    if(is.bipartite(x)){  #If bipartite, might need to adjust mode 1 count
      m1v<-get.network.attribute(x,"bipartite")  #How many mode 1 verts?
      set.network.attribute(x,"bipartite",m1v-sum(vid<=m1v))
    }
    x<-.Call(deleteVertices_R,x,vid)
    if(.validLHS(xn,parent.frame())){  #If x not anonymous, set in calling env 
      on.exit(eval.parent(call('<-',xn,x)))
    }
  }
  invisible(x)
}


# Retrieve a specified edge attribute from edge list el.  The attribute
# is returned as a list, unless unlist is TRUE. 
# if deleted.edges.omit is TRUE, then only attribute values on existing (non-null) edges will be returned.
# if na.omit is TRUE, than values corresponding to 'missing' edges (edges with attribute 'na' set to TRUE) should be ommited. (NULL edgs count as not-missing)
# If null.na is TRUE, then values corresponding to  edges for which the attribute name was never set will be set to NA.  Otherwise, they will be NULL, which means they will be included when unlist=TRUE 
#
get.edge.attribute<-function(el, attrname, unlist=TRUE,na.omit=FALSE,null.na=FALSE,deleted.edges.omit=FALSE){
  if (is.network(el)) el <- el$mel

  if (!is.list(el))
    stop("el must be a network object or a list.")

  if (!is.character(attrname))
    stop("attrname must be a character vector.")

  if (!is.logical(unlist) || !is.logical(na.omit) || !is.logical(null.na) ||
      !is.logical(deleted.edges.omit))
    stop("na.omit, null.na, deleted.edges.omit must be a logical vector.")

  edges <- .Call(getEdgeAttribute_R,el,attrname,na.omit,null.na,deleted.edges.omit)

  if(unlist)
    unlist(edges)
  else
    edges
}


# Retrieve a specified edge attribute from all edges in x.
#
get.edge.value<-function(x, attrname, unlist=TRUE, na.omit=FALSE, null.na=FALSE, deleted.edges.omit=FALSE){
  get.edge.attribute(x,attrname,unlist,na.omit,null.na,deleted.edges.omit)
}

# Retrieve the ID numbers for all edges incident on v, in network x.  
# Outgoing or incoming edges are specified by neighborhood, while na.omit 
# indicates whether or not missing edges should be omitted.  The return value
# is a vector of edge IDs.
#
get.edgeIDs<-function(x, v, alter=NULL, neighborhood=c("out","in","combined"), na.omit=TRUE){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("get.edgeIDs requires an argument of class network.")
  #Do some reality checking
  n<-network.size(x)
  if((v<1)||(v>n))
    return(numeric(0))
  if((!is.null(alter))&&((alter<1)||(alter>n)))
    return(numeric(0))
  #Retrieve the edges
  if(!is.directed(x))
    neighborhood="combined"       #If undirected, out==in==combined
  else
    neighborhood=match.arg(neighborhood)
  #Do the deed
  .Call(getEdgeIDs_R,x,v,alter,neighborhood,na.omit)
}


# Retrieve all edges incident on v, in network x.  Outgoing or incoming
# edges are specified by neighborhood, while na.omit indicates whether
# or not missing edges should be omitted.  The return value is a list of
# edges.
#
get.edges<-function(x, v, alter=NULL, neighborhood=c("out","in","combined"), na.omit=TRUE){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("get.edges requires an argument of class network.")
  #Do some reality checking
  n<-network.size(x)
  if((v<1)||(v>n))
    return(list())
  if((!is.null(alter))&&((alter<1)||(alter>n)))
    return(list())
  #Retrieve the edges
  if(!is.directed(x))
    neighborhood="combined"       #If undirected, out==in==combined
  else
    neighborhood=match.arg(neighborhood)
  #Do the deed
  .Call(getEdges_R,x,v,alter,neighborhood,na.omit)
}

# get the the edge ids associated with a set of dayds
# as defined by a vector of tails and heads vertex ids
get.dyads.eids<-function(x,tails,heads,neighborhood = c("out", "in", "combined")){
  if(length(tails)!=length(heads)){
    stop('heads and tails vectors must be the same length for get.dyads.eids')
  }
  if (any(heads>network.size(x) | heads<1) | any(tails>network.size(x) | tails<1)){
    stop('invalid vertex id in heads or tails vector')
  }
  neighborhood<-match.arg(neighborhood)
  if (!is.directed(x)){
    neighborhood = "combined"
  }
  lapply(seq_along(tails),function(e){
    eid<-get.edgeIDs(x,v = tails[e],alter=heads[e],neighborhood=neighborhood)
    if(length(eid)>1){
      eid<-NA
      warning('get.dyads.eids found multiple edge ids for dyad ',tails[e],',',heads[e],' NA will be returned')
    }
    eid
  })
}


# Given a network and a set of vertices, return the subgraph induced by those
# vertices (preserving all associated metadata); if given two such sets, 
# return the edge cut (along with the associated vertices and meta-data) as
# a bipartite network.
#
get.inducedSubgraph<-function(x, v, alters=NULL, eid=NULL){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("get.inducedSubgraph requires an argument of class network.")
  #Do some reality checking
  n<-network.size(x)
  
  # are we doing this via eids, or v and alters
  if (is.null(eid)){  # do checks for v and alters
    if((length(v)<1)||any(is.na(v))||any(v<1)||any(v>n))
      stop("Illegal vertex selection in get.inducedSubgraph")
    if(!is.null(alters)){
      if((length(alters)<1)||any(is.na(alters))||any(alters<1)||any(alters>n)|| any(alters%in%v))
        stop("Illegal vertex selection (alters) in get.inducedSubgraph")
    }
    if (!is.null(eid)){
      warning('eid argument to get.inducedSubgraph ignored when using v or alter argument')
    }
  } else { # do checks for eids
    if (!is.numeric(eid)){
      stop('eid must be a numeric vector of edge ids')
    }
    if (!missing(v)){
      warning('v argument to get.inducedSubgraph ignored when using eid argument')
    }
    if (!is.null(alters)){
      warning('alters argument to get.inducedSubgraph ignored when using eid argument')
    }
    # check that eids are valid
    if (any(!eid%in%valid.eids(x))){
      stop('eid argument contains non-valid edge ids')
    }
    
  }
  
  #Start by making a copy of our target network (yes, this can be wasteful)
  #TODO: in most cases, probably faster to create a new network and only copy over what is needed
  newNet<-network.copy(x)
  
  if (is.null(eid)){  # using v and alter
    #Now, strip out what is needed, and/or permute in the two-mode case
    if(is.null(alters)){                    #Simple case
      delete.vertices(newNet,(1:n)[-v])           #Get rid of everyone else
    }else{                                  #Really an edge cut, but w/vertices
      nv<-length(v)
      na<-length(alters)
      newids<-sort(c(v,alters))
      newv<-match(v,newids)
      newalt<-match(alters,newids)
      delete.vertices(newNet,(1:n)[-c(v,alters)])  #Get rid of everyone else
      permute.vertexIDs(newNet,c(newv,newalt))    #Put the new vertices first
      #Remove within-group edges
      for(i in 1:nv)
        for(j in (i:nv)[-1]){
          torem<-get.edgeIDs(newNet,i,alter=j,neighborhood="combined",na.omit=FALSE)
          if(length(torem)>0)
            delete.edges(newNet,torem)
        }
      for(i in (nv+1):(nv+na))
        for(j in (i:(nv+na))[-1]){
          torem<-get.edgeIDs(newNet,i,alter=j,neighborhood="combined",na.omit=FALSE)
          if(length(torem)>0)
            delete.edges(newNet,torem)
        }
      newNet%n%"bipartite"<-nv   #Set bipartite attribute
    }
  } else {  # using eids instead of v and alters
    # delete all the edges not in eid
    removeEid<-setdiff(valid.eids(newNet),eid)
    delete.edges(newNet,removeEid)
    # find the set of vertices incident on the remaining edges
    v<-unique(c(unlist(sapply(newNet$mel, "[[", "outl")),unlist(sapply(newNet$mel, "[[", "inl"))))
    removeV<-setdiff(seq_len(network.size(newNet)),v)
    delete.vertices(newNet,removeV)
  }
  #Return the updated object
  newNet
}


# Retrieve a specified network-level attribute from network x.  The attribute
# type depends on the underlying storage mode, and cannot be guaranteed.
#
get.network.attribute<-function(x,attrname,unlist=FALSE){
  x <- x$gal[[attrname]]
  if(unlist){unlist(x)}else{x}
}


# Retrieve the neighborhood of v in network x.  Depending on the value of 
# type, the neighborhood in question may be in, out, or the union of the two.
# The return value for the function is a vector containing vertex IDs.
#
get.neighborhood<-function(x, v, type=c("out","in","combined"), na.omit=TRUE){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("get.neighborhood requires an argument of class network.")
  #Do some reality checking
  n<-network.size(x)
  if((v<1)||(v>n))
    return(numeric(0))
  #Retrieve the edges
  if(!is.directed(x))
    type="combined"       #If undirected, out==in==combined
  else
    type=match.arg(type)
  #Do the deed
  .Call(getNeighborhood_R,x,v,type,na.omit)
}


# Retrieve a specified vertex attribute (indicated by attrname) from network x.
# Where na.omit==TRUE, values for missing vertices are removed; where
# null.na==TRUE, NULL values are converted to NAs.  The return value of this
# function is a list.
# 
get.vertex.attribute<-function(x,attrname,na.omit=FALSE,null.na=TRUE,
                               unlist=TRUE){
  #Check to see if there's anything to be done
  if(!is.network(x))
    stop("get.vertex.attribute requires an argument of class network.")
  if(network.size(x)==0){
    return(NULL)
  }
  #if(!(attrname %in% list.vertex.attributes(x))) 
  #  warning(paste('attribute', attrname,'is not specified for these vertices'))
  #Get the list of attribute values
  va<-lapply(x$val,"[[",attrname)
  #If needed, figure out who's missing
  if(na.omit)
    vna<-unlist(lapply(x$val,"[[","na"))
  else
    vna<-rep(FALSE,length(va))
  #Replace NULL values with NAs, if requested
  if(null.na)
    va[sapply(va,is.null)]<-NA
  #Return the result
  if (na.omit){
   x <- va[!vna]
  } else {
    x<-va
  }
  if(unlist){unlist(x)}else{x}
}


# Return TRUE iff network x has loops.
#
has.loops<-function(x){
  if(!is.network(x))
    stop("has.loops requires an argument of class network.")
  else
    get.network.attribute(x,"loops")
}


# Return TRUE iff (vi,vj) in network x.  Where na.omit==TRUE, edges flagged
# as missing are ignored.
#
is.adjacent<-function(x,vi,vj,na.omit=FALSE){
  if(!is.network(x))
    stop("is.adjacent requires an argument of class network.\n")
  if(length(vi)!=length(vj)){
    vi<-rep(vi,length=max(length(vi),length(vj)))
    vj<-rep(vj,length=max(length(vi),length(vj)))
  }
  #Do the deed
 .Call(isAdjacent_R,x,vi,vj,na.omit)
}


# Return TRUE iff network x is bipartite
#
is.bipartite<-function(x){
  if(!is.network(x))
    stop("is.bipartite requires an argument of class network.")
  else
    bip <- get.network.attribute(x,"bipartite")
  if(is.null(bip)){
   return(FALSE)
  } else if (is.logical(bip)){
   return(bip)  
  }else{
   return(bip>=0)
  }
}


# Return TRUE iff network x is directed.
#
is.directed<-function(x){
  if(!is.network(x))
    stop("is.directed requires an argument of class network.\n")
  else
    get.network.attribute(x,"directed")
}


# Return TRUE iff network x is hypergraphic.
#
is.hyper<-function(x){
  if(!is.network(x))
    stop("is.hyper requires an argument of class network.\n")
  else
    get.network.attribute(x,"hyper")
}


# Return TRUE iff network x is multiplex.
#
is.multiplex<-function(x){
  if(!is.network(x))
    stop("is.multiplex requires an argument of class network.\n")
  else
    get.network.attribute(x,"multiple")
}


# Return a network whose edges are the missing edges of x
#
is.na.network<-function(x){
  #Create an empty network with the same properties as x
  y<-network.initialize(network.size(x),directed=is.directed(x), hyper=is.hyper(x),loops=has.loops(x),multiple=is.multiplex(x), bipartite=x%n%"bipartite")
  #Add the missing edges of x to y
  y<-.Call(isNANetwork_R,x,y)
  #Return the updated network 
  y
}


# Return TRUE iff x is a network.
#
is.network<-function(x){
  inherits(x, "network")
}


# List attributes present on any edge
#
list.edge.attributes<-function(x){
  #First, check to see that this is a graph object
  if(!is.network(x))
    stop("list.edge.attributes requires an argument of class network.\n")
  # no edges in the network
  if (network.edgecount(x, na.omit=F) == 0) return(character(0))
  #Accumulate names
  allnam<-sapply(lapply(x$mel[!is.null(x$mel)],"[[","atl"),names)
  #Return the sorted, unique attribute names
  sort(unique(as.vector(unlist(allnam))))
}


# List network-level attributes
#
list.network.attributes<-function(x){
  #First, check to see that this is a graph object
  if(!is.network(x))
    stop("list.network.attributes requires an argument of class network.\n")
  #Return the attribute names
  sort(names(x$gal))
}


# List attributes present on any vertex
#
list.vertex.attributes<-function(x){
  #First, check to see that this is a graph object
  if(!is.network(x))
    stop("list.vertex.attributes requires an argument of class network.\n")
  if(network.size(x)==0){
    return(NULL)
  }
  #Accumulate names
  allnam<-unlist(sapply(x$val,names))
  #Return the sorted, unique attribute names
  sort(unique(as.vector(allnam)))
}


# Retrieve the number of free dyads (i.e., number of non-missing) of network x.
#
network.dyadcount<-function(x,na.omit=TRUE){
 if(!is.network(x)){
   stop("network.dyadcount requires an argument of class network.")
 }

 nodes <- network.size(x)
 if(is.directed(x)){
   if(is.bipartite(x)){ # directed bipartite
     nactor <- get.network.attribute(x,"bipartite")
     nevent <- nodes - nactor
     dyads <- nactor * nevent *2
   }else{ # directed unipartite
    dyads <- nodes * (nodes-1)
    if(has.loops(x)){
      # add in the diagonal
      dyads<-dyads+nodes
    }
   }
 }else{ # undirected
  if(is.bipartite(x)){ # undirected bipartite
   nactor <- get.network.attribute(x,"bipartite")
   nevent <- nodes - nactor
   dyads <- nactor * nevent
  }else{ # undirected unipartite
   dyads <- nodes * (nodes-1)/2
   if(has.loops(x)){
     # add in the diagonal
     dyads<-dyads+nodes
   }
  }
 }
 if(na.omit){
#
#  Adjust for missing
#
  design <- get.network.attribute(x,"design")
  if(!is.null(design)){
   dyads <- dyads - network.edgecount(design)
  }else{
   design <- get.network.attribute(x,"mClist.design")
   if(!is.null(design)){
    dyads <- dyads - design$nedges
   }else{
    dyads <- dyads - network.naedgecount(x)
   }
  }
 }
 dyads
}


#Retrieve the number of edges in network x.
#
network.edgecount<-function(x,na.omit=TRUE){
  #First, check to see that this is a graph object
  if(!is.network(x))
    stop("network.edgecount requires an argument of class network.\n")
  #Return the edge count
  .Call(networkEdgecount_R,x,na.omit)
}


#Retrieve the number of missing edges in network x
#
network.naedgecount<-function(x){
  na<-get.edge.attribute(x$mel,"na")
  if(is.null(na))
    0
  else
    sum(na)
}


# Retrieve the size (i.e., number of vertices) of network x.
#
network.size<-function(x){
  if(!is.network(x))
    stop("network.size requires an argument of class network.\n")
  else
    get.network.attribute(x,"n")
}


# Retrieve the vertex names of network x (if present).
#
network.vertex.names<-function(x){
  if(!is.network(x)){
    stop("network.vertex.names requires an argument of class network.")
  }else{
    if(network.size(x)==0)
      return(NULL)
    vnames <- get.vertex.attribute(x,"vertex.names")
    if(is.null(vnames)  | all(is.na(vnames)) ){
      paste(1:network.size(x))
    }else{
      vnames
    }
  }
}


# Set the vertex names of network x
#
"network.vertex.names<-"<-function(x,value){
  set.vertex.attribute(x,attrname="vertex.names",value=value)
}


# Permute the internal IDs (ordering) of the vertex set
permute.vertexIDs<-function(x,vids){
  #First, check to see that this is a graph object
  if(!is.network(x))
    stop("permute.vertexIDs requires an argument of class network.\n")
  #Sanity check: is this a permutation vector?
  n<-network.size(x)
  if((length(unique(vids))!=n)||(range(vids)!=c(1,n)))
    stop("Invalid permutation vector in permute.vertexIDs.")
  if(is.bipartite(x)){  #If bipartite, enforce partitioning
    bpc<-get.network.attribute(x,"bipartite")
    if(any(vids[0:bpc]>bpc)||(vids[(bpc+1):n]<=bpc))
      warning("Performing a cross-mode permutation in permute.vertexIDs.  I hope you know what you're doing....")
  }
  #Return the permuted graph
  xn<-substitute(x)
  x<-.Call(permuteVertexIDs_R,x,vids)
  if(.validLHS(xn,parent.frame())){  #If x not anonymous, set in calling env 
    on.exit(eval.parent(call('<-',xn,x)))
  }
  invisible(x)
}


# Set an edge attribute for network x.
#
# set.edge.attribute<-function(x,attrname,value,e=seq_along(x$mel)){
#   #Check to be sure we were called with a network
#   if(!is.network(x))
#     stop("set.edge.attribute requires an argument of class network.")
#   #Make sure that value is appropriate, coercing if needed
#   if(!is.list(value)){
#     if(!is.vector(value))
#       stop("Inappropriate edge value given in set.edge.attribute.\n")
#     else
#       value<-as.list(rep(value,length=length(e)))
#   }else
#     if(length(value)!=length(e))
#       value<-rep(value,length=length(e))
#   xn<-deparse(substitute(x))
#   ev<-parent.frame()
#   if(length(e)>0){
#     if((min(e)<1)|(max(e)>length(x$mel)))
#       stop("Illegal edge in set.edge.attribute.\n")
#     #Do the deed
#     x<-.Call("setEdgeAttribute_R",x,attrname,value,e, PACKAGE="network")
#     if(exists(xn,envir=ev))          #If x not anonymous, set in calling env
#       on.exit(assign(xn,x,pos=ev))
#     invisible(x)
#   }else
#     invisible(x)
# }

set.edge.attribute<-function(x,attrname,value,e=seq_along(x$mel)){
  #Check to be sure we were called with a network
  if(!is.network(x)){
    stop("set.edge.attribute requires an argument of class network.")
  }
  # determine if we have to do anything at all
  if(length(e)>0){
    if((min(e)<1)|(max(e)>length(x$mel))){
      stop("Illegal edge in set.edge.attribute.\n")
    }
    xn<-substitute(x)
    # determine if we will be setting single or multiple values
    if(length(attrname)==1){
      #Make sure that value is appropriate, coercing if needed
      if(!is.list(value)){
        if(!is.vector(value)){
          stop("Inappropriate edge value given in set.edge.attribute.\n")
        } else {
          value<-as.list(rep(value,length=length(e)))
        }
      } else {
        if(length(value)!=length(e)) {
          value<-rep(value,length=length(e))
        }
      }
      #Do the deed, call the set single value version
      x<-.Call(setEdgeAttribute_R,x,attrname,value,e)
    } else { # we will be setting multiple values
      if (length(attrname)!=length(value)){
        stop("the 'value' attribute must have an element corresponding to each attribute name in 'attrname' in set.edge.attribute")
      }
      #Make sure that value is appropriate, coercing if needed
      if(!is.list(value)){
        if(!is.vector(value)){
          stop("Inappropriate edge value given in set.edge.attribute.\n")
        } else { # value must be a vector
          # replicate each element of value e times if needed
          value<-lapply(1:length(value),function(n){
            if (length(value[n])<length(e)){
              return(as.list(rep(value[n],length=length(e))))
            } else {
              return(as.list(value[n]))
            }
          })
        }
      } else {
        # replicate each element of value e times if needed
        value<-lapply(1:length(value),function(n){
          if (length(value[[n]])<length(e)){
            return(as.list(rep(value[[n]],length=length(e))))
          } else {
            return(as.list(value[[n]]))
          }
        })
        
      }
      #Do the deed, call the set multiple version
      x<-.Call(setEdgeAttributes_R,x,attrname,value,e)
    }
    if(.validLHS(xn,parent.frame())){  #If x not anonymous, set in calling env 
      on.exit(eval.parent(call('<-',xn,x)))
    }
  }
  invisible(x)
}


# Set an edge value for network x.
#
set.edge.value<-function(x,attrname,value,e=seq_along(x$mel)){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("set.edge.value requires an argument of class network.\n")
  #Check to ensure that this is not a hypergraph
  if(is.hyper(x))
    stop("Hypergraphs not currently supported in set.edge.value.\n")
  # Check edges
  if (length(e)==0) return(invisible(x))
  if((min(e)<1)|(max(e)>length(x$mel)))
    stop("Illegal edge in set.edge.value.\n")
  #Make sure that value is appropriate, coercing if needed
  n<-network.size(x)
  if(!is.matrix(value)){
    if(is.vector(value))
      value<-matrix(rep(value,length=n*n),n,n)
    else
      value<-matrix(value,n,n)
  } else if (min(dim(value)) < n) {
    stop("set.edge.value requires a matrix whose dimension is equal to or larger than the network size")
  }
  #Do the deed
  xn<-substitute(x)
  x<-.Call(setEdgeValue_R,x,attrname,value,e)
  if(.validLHS(xn,parent.frame())){  #If x not anonymous, set in calling env 
    on.exit(eval.parent(call('<-',xn,x)))
  }
  invisible(x)
}


# Set a network-level attribute for network x.
#
set.network.attribute<-function(x,attrname,value){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("set.network.attribute requires an argument of class network.")
  #Make sure the values are consistent
  if(length(attrname)==1){
    value<-list(value)
  }else{
    if(is.list(value)){
      value<-rep(value,length=length(attrname))
    }else if(is.vector(value)){
      value<-as.list(rep(value,length=length(attrname)))
    }else
      stop("Non-replicable value with multiple attribute names in set.network.attribute.\n")
  }
  #Do the deed
  xn<-substitute(x)
  x<-.Call(setNetworkAttribute_R,x,attrname,value)
  if(.validLHS(xn,parent.frame())){  #If x not anonymous, set in calling env 
    on.exit(eval.parent(call('<-',xn,x)))
  }
  invisible(x)
}


# Set a vertex attribute for network x.
# This version has been removed so we can test one that can set multiple values at once
# set.vertex.attribute<-function(x,attrname,value,v=seq_len(network.size(x))){
#   #Check to be sure we were called with a network
#   if(!is.network(x))
#     stop("set.vertex.attribute requires an argument of class network.")
#   #Perform some sanity checks
#   if(any((v>network.size(x))|(v<1)))
#     stop("Vertex ID does not correspond to actual vertex in set.vertex.attribute.\n")
#   #Make sure that value is appropriate, coercing if needed
#   if(!is.list(value)){
#     if(!is.vector(value))
#       stop("Inappropriate value given in set.vertex.attribute.\n")
#     else
#       value<-as.list(rep(value,length=length(v)))
#   }else
#     if(length(value)!=length(v))
#       value<-rep(value,length=length(v))
#   #Do the deed
#   xn<-deparse(substitute(x))
#   ev<-parent.frame()
#   x<-.Call("setVertexAttribute_R",x,attrname,value,v, PACKAGE="network")
#   if(exists(xn,envir=ev))          #If x not anonymous, set in calling env
#     on.exit(assign(xn,x,pos=ev))
#   invisible(x)
# }

# valid.eids  returns a list of non-null edge ids for a given network
valid.eids <-function(x){
  # maybe should omit class test for speed?
  if (!is.network(x)){
    stop("cannot determine non-null edge ids because argument x is not a network object")
  }
  # get the ids of all the non-null elements on the edgelist of x
  return(which(!sapply(x$mel,is.null)))
}

set.vertex.attribute<-function(x,attrname,value,v=seq_len(network.size(x))){
  #Check to be sure we were called with a network
  if(!is.network(x))
    stop("set.vertex.attribute requires an argument of class network.")
  #Perform some sanity checks
  if(any((v>network.size(x))|(v<1)))
    stop("Vertex ID does not correspond to actual vertex in set.vertex.attribute.\n")
  
  xn<-substitute(x)
  
  #Make sure that value is appropriate, coercing if needed
  if (length(attrname)==1){ # if we are only setting a single attribute use old version
    if(!is.list(value)){
      if(!is.vector(value)){
        stop("Inappropriate value given in set.vertex.attribute.\n")
      } else {
        value<-as.list(rep(value,length=length(v)))
      }
    } else {
      if(length(value)!=length(v)){
        value<-rep(value,length=length(v))
      }
    }
    # call older singular value version
    x<-.Call(setVertexAttribute_R,x,attrname,value,v)
  } else { # setting multiple values
    if (length(value)!=length(attrname)){
      stop("the 'value' attribute must have an element corresponding to each attribute name in 'attrnames' in set.vertex.attribute")
    }
    if(!is.list(value)){
      if(!is.vector(value)){
        stop("Inappropriate value given in set.vertex.attribute.\n")
      } else { # value is a vector
    
        # replicate each element of value v times if needed
        value<-lapply(1:length(value),function(n){
                  if (length(value[n])<length(v)){
                    return(as.list(rep(value[n],length=length(v))))
                  } else {
                    return(as.list(value[n]))
                  }
              })
      }
    } else {  # value is a list
      # replicate each element of value v times if needed
      value<-lapply(1:length(value),function(n){
        if (length(value[[n]])<length(v)){
          return(as.list(rep(value[[n]],length=length(v))))
        } else {
          return(as.list(value[[n]]))
        }
      })
    }
    # call multiple value version
    x<-.Call(setVertexAttributes_R,x,attrname,value,v)
  } # end setting multiple values
  #Do the deed
  
  if(.validLHS(xn,parent.frame())){  #If x not anonymous, set in calling env 
    on.exit(eval.parent(call('<-',xn,x)))
  }
  invisible(x)
}

# valid.eids  returns a list of non-null edge ids for a given network
valid.eids <-function(x){
  # maybe should omit class test for speed?
  if (!is.network(x)){
    stop("cannot determine non-null edge ids because argument x is not a network object")
  }
  # get the ids of all the non-null elements on the edgelist of x
  return(which(!sapply(x$mel,is.null)))
}

