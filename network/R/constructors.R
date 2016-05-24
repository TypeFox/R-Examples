######################################################################
#
# constructors.R
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
# This file contains various routines for the construction of network
# and edge objects.
#
# Contents:
#
#   network
#   network.adjacency
#   network.copy
#   network.edgelist
#   network.incidence
#   network.initialize
#
######################################################################


# Basic network constructor.  Converts a single matrix to a network class
# object.  The matrix must be in one of three formats:  adjacency,
# incidence, or edgelist.
#
# MSH added bipartite
#

network<-function(x, vertex.attr=NULL, vertex.attrnames=NULL,
                directed=TRUE, hyper=FALSE, loops=FALSE,
                multiple=FALSE, bipartite=FALSE, ...)
{
  #Initialize the network object
  g<-as.network(x,directed=directed,hyper=hyper,loops=loops,
              multiple=multiple,bipartite=bipartite,...)
  #Add vertex attributes, if needed
  if(!is.null(vertex.attr)){
    #Create vertex attribute names, if needed
    if(is.null(vertex.attrnames)){
      if(!is.null(names(vertex.attr)))
        vertex.attrnames<-names(vertex.attr)
      else{
        vertex.attrnames<-1:length(vertex.attr)
	warning("Vertex attribute names not given; making some up.")
      }
    }
    #Add the attributes
    for(i in 1:length(vertex.attr))
      g<-set.vertex.attribute(g,vertex.attrnames[[i]],vertex.attr[[i]])
  }
# xnames <- get.vertex.attribute(g,"vertex.names")
# if(!is.null(xnames) & any(!is.na(xnames))){ g <- xnames }
  #Return the result
  g  
}

# Construct a network's edge set, using an a bipartite adjacency matrix as input.
#
network.bipartite<-function(x, g, ignore.eval=TRUE, names.eval=NULL, ...){
  #Set things up to edit g in place
  gn<-substitute(g)
  #Build head/tail lists; note that these cannot be hypergraphic or
  #multiplex, since our data is drawn from an adjacency matrix
  nactors <- dim(x)[1]
  nevents <- dim(x)[2]
  n <- nactors + nevents
  #Add names if available
  if(!is.null(colnames(x)) & !is.null(rownames(x))){
    g <- set.vertex.attribute(g,"vertex.names",c(rownames(x),colnames(x)))
  }
  # convert x into a matrix
  x<-as.matrix(x)
  
  X <- matrix(0,ncol=n,nrow=n)
# diag(X) <- 0
  X[1:nactors, nactors+(1:nevents)] <- x
  X[nactors+(1:nevents), 1:nactors] <- t(x)
  X[row(X)<col(X)]<-0            #Clear above-diagonal entries.
  x <- X
  missing <- is.na(x)
  x[missing] <- 1
#
  x<-as.vector(x)
  n<-network.size(g)
  e<-(0:(n*n-1))[x!=0] 
  if(ignore.eval){
    ev<-as.list(as.logical(missing[x!=0]))
    en<-replicate(length(ev),list("na"))
  }else{
    xv<-x
    ev<-apply(cbind(as.list(as.logical(missing[x!=0])),as.list(xv[x!=0])),1, as.list)
    en<-replicate(length(ev),list(list("na",names.eval)))
  }
  if(sum(x!=0)>0)
    add.edges(g, as.list(1+e%%n), as.list(1+e%/%n),
              names.eval=en, vals.eval=ev, ...)
  #Patch up g on exit for in-place modification
  if(.validLHS(gn,parent.frame())){
    on.exit(eval.parent(call('<-',gn,g)))
  }
  invisible(g)
}


# Construct a network's edge set, using an adjacency matrix as input.
#
network.adjacency<-function(x, g, ignore.eval=TRUE, names.eval=NULL, ...){
  # check that dimension of g is appropriate for x
  if (nrow(x)!=ncol(x)){
    stop('the network.adjacency constructor expects its matrix argument to be square (same number of rows and columns)')
  }
  if (network.size(g) != nrow(x)){
    stop('the network.adjacency constructor requires that the size of its network argument (',network.size(g),') matches the dimensions of the matrix argument (',nrow(x),' by ',ncol(x),')')
  }
  
  #Set things up to edit g in place
  gn<-substitute(g)
  #Build head/tail lists; note that these cannot be hypergraphic or
  #multiplex, since our data is drawn from an adjacency matrix
  if(!is.directed(g)){
    missingE <- is.na(x) | is.na(t(x))
    x[missingE] <- 1
    #Be sure to pick up nonzero entries for which x[i,j]=-x[j,i].
    x[x==-t(x)]<-abs(x)[x==-t(x)]  
    x<-(x+t(x))/2                  #Symmetrize matrix.
    x[row(x)<col(x)]<-0            #Clear above-diagonal entries.
  }else{
    missingE <- is.na(x)
    x[missingE] <- 1
  }
  
  # if the na.rm value is specified and TRUE, don't include those missing edges after all
  # some ugliness to pull names from ...
  dotNames<-as.list(substitute(list(...)))[-1L]
  if('na.rm'%in%dotNames){
    na.rm<-list(...)[[match('na.rm',dotNames)]]
    if (na.rm){
      x[missingE]<-0
    }
  }
  
  if(!has.loops(g)){ # if it doesn't have loops, replace the diagonal
    diag(x)<-0
  }
  x<-as.vector(x)
  n<-network.size(g)
  e<-(0:(n*n-1))[x!=0] 
  if(ignore.eval){
    ev<-as.list(as.logical(missingE[x!=0]))
    en<-replicate(length(ev),list("na"))
  }else{
    xv<-x
    xv[missingE]<-NA
    ev<-apply(cbind(as.list(as.logical(missingE[x!=0])),as.list(xv[x!=0])),1, as.list)
    en<-replicate(length(ev),list(c("na",names.eval)))
  }
  # Add names if available
  if(!is.null(colnames(x))){
   g <- set.vertex.attribute(g,"vertex.names", colnames(x))
  }else{
    if(!is.null(rownames(x))){
      g <- set.vertex.attribute(g,"vertex.names", rownames(x))
    }
  }
  if(sum(x!=0)>0)
    add.edges(g, as.list(1+e%%n), as.list(1+e%/%n),
              names.eval=en, vals.eval=ev, ...)
  #Patch up g on exit for in-place modification
  if(.validLHS(gn,parent.frame())){
    on.exit(eval.parent(call('<-',gn,g)))
  }
  invisible(g)
}


# Construct and a return a network object which is a copy of x
#
network.copy<-function(x){
  #Verify that this is a network object
  if(!is.network(x))
    stop("network.copy requires an argument of class network.\n")
  #Duplicate and return
  y<-.Call(copyNetwork_R,x)
  y
}


# Construct a network's edge set, using an edgelist matrix as input.
#
network.edgelist<-function(x, g, ignore.eval=TRUE, names.eval=NULL, ...){
  #Set things up to edit g in place
  gn<-substitute(g)
  l<-dim(x)[2]
  #Traverse the edgelist matrix, adding edges as we go.
  if((l>2)&&(!ignore.eval)){		#Use values if present...
    #if names not given, try to use the names from data frame
    if (is.null(names.eval)){
      names.eval<-names(x)[3:l]
    }
    #if it is still null, its going to crash, so throw an informative error
    if (is.null(names.eval)){
      stop("unable to add attribute values to edges because names are not provided for each attribute (names.eval=NULL)")
    }
    edge.check<-list(...)$edge.check 
    eattrnames <-lapply(seq_len(NROW(x)),function(r){as.list(names.eval)})
   # eattrvals <-apply(x[,3:l,drop=FALSE]
    eattrvals <-lapply(seq_len(NROW(x)),function(r){as.list(x[r,3:l,drop=FALSE])})
    g<-add.edges(g,as.list(x[,1]),as.list(x[,2]),eattrnames,eattrvals,edge.check=edge.check)
  }else{				#...otherwise, don't.
    edge.check<-list(...)$edge.check      
    g<-add.edges(g,as.list(x[,1]),as.list(x[,2]),edge.check=edge.check)
  }
  #Patch up g on exit for in-place modification
  if(.validLHS(gn,parent.frame())){
    on.exit(eval.parent(call('<-',gn,g)))
  }
  invisible(g)
}


# Construct a network's edge set, using an incidence matrix as input.
#
network.incidence<-function(x, g, ignore.eval=TRUE, names.eval=NULL, ...){
  #Set things up to edit g in place
  gn<-substitute(g)
  n<-network.size(g)
  edge.check<-list(...)$edge.check      
  #Traverse the incidence matrix, adding edges as we go.
  for(i in 1:dim(x)[2]){
    #Construct the head and tail sets
    if(is.directed(g)){
      if(any(is.na(x[,i])))
        stop("Missing data not allowed for directed incidence matrices.\n")
      head<-(1:n)[x[,i]>0]
      tail<-(1:n)[x[,i]<0]
      missing<-FALSE
    }else{
      missing<-any(is.na(x[,i]))
      x[,i][is.na(x[,i])]<-1
      head<-(1:n)[x[,i]!=0]
      if(is.hyper(g))
        tail<-head
      else{                 #If dyadic, use only the first two nonzero entries
        tail<-head[1]
        head<-head[2]
      }
    }
    if(length(head)*length(tail)==0)
      stop("Supplied incidence matrix has empty head/tail lists. (Did you get the directedness right?)")
    #Get edge values, if needed
    if(ignore.eval){
      en<-"na"
      ev<-missing
    }else{
      if(!is.directed(g))
        ev<-list(missing,x[x[,i]!=0,i][1])
      else
        ev<-list(missing,abs(x[x[,i]!=0,i][1]))
      if(is.null(names.eval))
        en<-list("na",NULL)
      else
        en<-list("na",names.eval)
    }
    #Add the edge to the graph
    g<-add.edge(g,tail,head,names.eval=en,vals.eval=ev,edge.check=edge.check)
  }
  #Patch up g on exit for in-place modification
  if(.validLHS(gn,parent.frame())){
    on.exit(eval.parent(call('<-',gn,g)))
  }
  invisible(g)
}

# Initialize a new network object.
# MSH added bipartite
#
network.initialize<-function(n,directed=TRUE,hyper=FALSE,loops=FALSE,multiple=FALSE,bipartite=FALSE){
  #If we have a negative number of vertices, we have a problem...
  n<-round(n)
  if(n<0)
    stop("Network objects cannot be of negative order.")
  #Create the base-level lists
  g<-list()
  g$mel<-list()
  g$gal<-list()
  #Create the required network attributes
  g$gal$n<-n
  g$gal$mnext<-1
  g$gal$directed<-directed
  g$gal$hyper<-hyper
  g$gal$loops<-loops
  g$gal$multiple<-multiple
  g$gal$bipartite<-bipartite
  #Populate the vertex attribute lists, endpoint lists, etc.
  if(n>0){
    g$val<-replicate(n,list())
    g$iel<-replicate(n,vector(mode="integer"))
    g$oel<-replicate(n,vector(mode="integer"))
  }else{
    g$val<-vector(length=0,mode="list")
    g$iel<-vector(length=0,mode="list")
    g$oel<-vector(length=0,mode="list")
  }
  #Set the class
  class(g)<-"network"
  #Set the required vertex attribute
  if(n>0)
    g<-set.vertex.attribute(g,"na",rep(FALSE,n),1:n)
  #Create default vertex names
  if(n>0)
    network.vertex.names(g)<-1:n
  #Return
  g
}
