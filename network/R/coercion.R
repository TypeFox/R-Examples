######################################################################
#
# coercion.R
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
# This file contains various routines for coercion to/from network
# class objects.
#
# Contents:
#
#   as.matrix.network
#   as.matrix.network.adjacency
#   as.matrix.network.edgelist
#   as.matrix.network.incidence
#   as.network
#   as.network.default
#   as.network.network
#   as.network.matrix
#   as.sociomatrix
#
######################################################################


# Method for general coercion of network class objects into matrices.
# Matrix type is indicated by the eponymous argument; note that some
# types may not be supported for certain networks.  Where
# attrname!=NULL, an edge attribute of name attrname is used to supply
# edge values.  Otherwise, edges are assumed to be unvalued.
#
as.matrix.network<-function(x,matrix.type=NULL,attrname=NULL,...){
  #Get the matrix type
  if(is.null(matrix.type))
    matrix.type<-"adjacency"
  else
    matrix.type<-match.arg(matrix.type,c("adjacency","incidence","edgelist"))
  #Dispatch as needed
  switch(matrix.type,
    adjacency=as.matrix.network.adjacency(x=x,attrname=attrname,...),
    incidence=as.matrix.network.incidence(x=x,attrname=attrname,...),
    edgelist=as.matrix.network.edgelist(x=x,attrname=attrname,...)
  )
}


# Coerce a network object to an adjacency matrix (where possible).  If
# provided, attrname is used to identify an attribute to use for edge
# values.
#
as.matrix.network.adjacency<-function(x,attrname=NULL,expand.bipartite=FALSE,...){
  #Check to make sure this is a supported network type
  if(is.hyper(x))
    stop("Hypergraphs not currently supported in as.matrix.network.adjacency.  Exiting.\n")
  if(is.multiplex(x))
    stop("Multigraphs not currently supported in as.matrix.network.adjacency.  Exiting.\n")
  #Generate the adjacency matrix 
  m<-matrix(0,nrow=network.size(x),ncol=network.size(x))
  if(network.size(x)==0)
    return(m)
  tl<-unlist(sapply(x$mel,"[[","outl")) #Can unlist b/c no hyperedges
  hl<-unlist(sapply(x$mel,"[[","inl"))
  nal<-as.logical(get.edge.attribute(x$mel,"na",unlist=TRUE))
  if(!is.null(attrname)){
    val<-unlist(get.edge.attribute(x$mel,attrname,unlist=FALSE))
    if(is.null(val)){
     warning(paste("There is no edge attribute named", attrname))
     val<-rep(1,length(tl))
    }
  }else{
    val<-rep(1,length(tl))
  }
  if(length(hl[!nal])>0){
    m[tl[!nal]+(hl[!nal]-1)*network.size(x)]<-val[!nal]
  }
  if(length(hl[ nal])>0){
   m[tl[ nal]+(hl[ nal]-1)*network.size(x)]<-NA
  }
  #If undirected, symmetrize
  if(!is.directed(x)){
# changed by MSH to allow non binary values
#   m<-pmax(m,t(m))
    sel<-m
    sel[is.na(m)]<-1
    m[sel==0] <- t(m)[sel==0]
  }
  #Set row/colnames to vertex names
  xnames <- network.vertex.names(x)
  dimnames(m) <- list(xnames, xnames)
  #If bipartite and !expand.bipartite, return in two-mode form
  if(is.bipartite(x)&(!expand.bipartite)){
    nactors <- get.network.attribute(x, "bipartite")
    nevents <- network.size(x) - nactors
    m <- m[0:nactors, nactors+(1:nevents)]
  }
  #Return the result
  m
}


# Coerce a network object to an edgelist matrix.  If provided, attrname is 
# used to identify an attribute to use for edge values.  Setting as.sna.edgelist
# results in output in the sna edgelist format (including missing edge handling)
# and is used by the sna package for coercion.
#
as.matrix.network.edgelist<-function(x,attrname=NULL,as.sna.edgelist=FALSE,...){
  #Check to make sure this is a supported network type
  if(is.hyper(x))
    stop("Hypergraphs not currently supported in as.matrix.network.edgelist.  Exiting.\n")
  #Find the missing edges
  nal<-as.logical(get.edge.attribute(x$mel,"na"))
  #Generate the edgelist matrix
  m<-cbind(unlist(sapply(x$mel,"[[","outl")), unlist(sapply(x$mel,"[[","inl")))
  #Add edge values, if needed
  if(!is.null(attrname))
    m<-cbind(m,get.edge.attribute(x$mel,attrname,na.omit=FALSE,null.na=TRUE,deleted.edges.omit=TRUE))
  else if(as.sna.edgelist)
    m<-cbind(m,rep(1,NROW(m)))
  #Set additional attributes and return the result
  if(as.sna.edgelist && nrow(m) > 0) # check that there are actually edges
    m[nal,3]<-NA
  else
    m<-m[!nal,,drop=FALSE]
  if(length(m)==0)
    m<-matrix(numeric(0),ncol=2+as.sna.edgelist)
  else if((!is.directed(x))&&as.sna.edgelist){    #sna uses directed form
    m<-rbind(m,m[m[,2]!=m[,1],c(2:1,3)])
  }
  attr(m,"n")<-network.size(x)
  attr(m,"vnames")<-network.vertex.names(x)
  if(is.bipartite(x))
    attr(m,"bipartite")<-x%n%"bipartite"
  m
}


# Coerce a network object to an incidence matrix (where possible).  If
# provided, attrname is used to identify an attribute to use for edge
# values.
#
as.matrix.network.incidence<-function(x,attrname=NULL,...){
  #Perform preprocessing
  n<-network.size(x)
  nulledge<-sapply(x$mel,is.null)
  inl<-lapply(x$mel,"[[","inl")[!nulledge]
  outl<-lapply(x$mel,"[[","outl")[!nulledge]
  if(!is.null(attrname))
    evals<-unlist(get.edge.attribute(x$mel,attrname))[!nulledge]
  else
    evals<-rep(1,length(x$mel))[!nulledge]
  ena<-as.logical(get.edge.attribute(x$mel,"na"))[!nulledge]
  #If called with an empty graph, return a degenerate matrix
  if(length(ena)==0)
    return(matrix(numeric(0),nrow=n))
  #Generate the incidence matrix
  dir<-is.directed(x)
  f<-function(a,m,k){y<-rep(0,m); y[a]<-k; y}
  im<-sapply(inl,f,n,1)+sapply(outl,f,n,ifelse(dir,-1,1))
  if(!dir)
    im<-pmin(im,1)
  im<-sweep(im,2,evals,"*")              #Fill in edge values
  im[(sapply(ena,rep,n)*(im!=0))>0]<-NA      #Add NAs, if needed
  #Return the result
  im
}


as.network<-function(x,...)
  UseMethod("as.network")


as.network.default<-function(x,...)
  as.network.matrix(x,...)


as.network.network<-function(x,...)
  x


#
# MSH modified for bipartite
#
as.network.matrix<-function(x, matrix.type=NULL,
        directed=TRUE, hyper=FALSE, loops=FALSE, multiple=FALSE,
        bipartite=FALSE,
        ignore.eval=TRUE, names.eval=NULL, na.rm=FALSE, edge.check=FALSE, ...){
  #Before doing anything else, pull any attributes from the matrix that we
  #might need....
  nattr<-attr(x,"n")             #Currently, only using sna edgelist attributes
  battr<-attr(x,"bipartite")
  vattr<-attr(x,"vnames")
  #Convert logicals to numeric form
  if(is.logical(x)){x <- 1*x}
  #Get the matrix type
  if(is.null(matrix.type))
    matrix.type<-which.matrix.type(x)
  else
    matrix.type<-match.arg(matrix.type,c("adjacency","incidence","edgelist",
                                         "bipartite"))
  if(is.logical(bipartite)&&bipartite)
    matrix.type<-"bipartite"
  #Patch adj->bipartite case
  if((bipartite>0)&&(matrix.type=="adjacency")&&(NROW(x)==bipartite))  
    matrix.type<-"bipartite"
  
  # Add names if available
  unames <- NULL
  if(matrix.type=="edgelist"){
    if(dim(x)[2]>2)
      vals<-x[,-(1:2),drop=FALSE]
    else
      vals<-NULL
    if(is.character(x<-as.matrix(x[,1:2,drop=FALSE]))){
      unames <- sort(unique(as.vector(x)))
      x <- cbind(match(x[,1],unames),match(x[,2],unames))
    }
    if(!is.null(vals)){
      x<-cbind(x,vals)
      
      if (is.null(colnames(vals))){
        colnames(x)<-NULL  #R creates these, and they are annoying later
      } else {
        # leave colnames for vals intact so they can be used for edge attributes
        colnames(x)<-c(NA,NA,colnames(vals))
      }
    }
  }
  if(matrix.type=="adjacency" && !is.null(colnames(x))){
    unames <- colnames(x)
  }
  if(matrix.type=="bipartite"){
   directed <- FALSE
   bipartite <- dim(x)[1]
   unames <- 1:sum(dim(x))
   if(!is.null(rownames(x))){
     unames[1:(dim(x)[1])] <- rownames(x)
   }
   if(!is.null(colnames(x))){
     unames[(dim(x)[1])+(1:(dim(x)[2]))] <- colnames(x)
   }
  }
  if(!is.null(vattr))                        #If given names, use 'em
    unames<-vattr
  #Initialize the network object
  n<-switch(matrix.type,	#Extract n based on matrix type
    adjacency=dim(x)[1],
    incidence=dim(x)[1],
    bipartite=sum(dim(x)),
    edgelist=max(x[,1:2]),
  )
  if(is.numeric(nattr))                      #If given n, use it
    n<-nattr
  if(is.numeric(battr))                      #If given bipartite info, use it
    bipartite<-battr
  
  # if we are going to build an adjacency matrix and it doesn't match the nattr, give an error, because otherwise will crash
  # this may happen if a square edgelist with attribute information is passed in
  if (is.numeric(nattr) & matrix.type=='adjacency'){
    if (nattr != ncol(x)){
      stop('the dimensions of the matrix argument (',nrow(x),' by ', ncol(x),') do not match the network size indicated by the attached n attribute (',nattr,'), perhaps matrix.type argument is not correct')
    }
  }
  
  g<-network.initialize(n,directed=directed, hyper=hyper, loops=loops, multiple=multiple,bipartite=bipartite)
  #Call the specific coercion routine, depending on matrix type
  g<-switch(matrix.type,
    adjacency=network.adjacency(x,g,
     ignore.eval,names.eval,na.rm,edge.check),
    incidence=network.incidence(x,g,
     ignore.eval,names.eval,na.rm,edge.check),
    bipartite=network.bipartite(x,g,
     ignore.eval,names.eval,na.rm,edge.check),
    edgelist=network.edgelist(x,g, 
     ignore.eval,names.eval,na.rm,edge.check)
  )

  if(!is.null(unames)){
   g <- set.vertex.attribute(g,"vertex.names", unames)
  }
  #Return the result
  g
}


#Force the input into sociomatrix form.  This is a shortcut to 
#as.matrix.network.adjacency, which ensures that a raw matrix is
#passed through as-is.
as.sociomatrix<-function(x, attrname=NULL, simplify=TRUE, expand.bipartite=FALSE, ...){
  if(is.network(x)){ #If network, coerce to adjacency matrix
    g<-as.matrix.network.adjacency(x,attrname=attrname, expand.bipartite=expand.bipartite,...)
  }else if(is.matrix(x)||is.array(x)){ #If an array/matrix, use as-is
    g<-x
  }else if(is.list(x)){  #If a list, recurse on list elements
    g<-lapply(x,as.sociomatrix,attrname=attrname,simplify=simplify, expand.bipartite=expand.bipartite,...)
  }else{
    stop("as.sociomatrix input must be an adjacency matrix/array, network, or list.")
  }
  #Convert into the appropriate return format
  if(is.list(g)){   #Collapse if needed
    if(length(g)==1){
      g<-g[[1]]
      if((!simplify)&&(length(dim(g))==3)){  #Coerce to a list of matrices?
        out<-list()
        for(i in 1:dim(g)[1])
          out[[i]]<-g[i,,]
      }else{
        out<-g
      }
    }else{
      #Coerce to array form?
      if(simplify){
        dims<-sapply(g,dim)
        if(is.list(dims)){      #Dims must not be of equal length
          mats<-sapply(dims,length)
          mats[mats==1]<-0
          mats[mats==2]<-1
          mats[mats==3]<-sapply(dims[mats==3],"[[",1)
          mats<-cumsum(mats)
          dims<-sapply(dims,"[",2)
        }else{                  #Dims are of equal length
          if(NROW(dims)==3)      #Determine number of matrices per entry
            mats<-cumsum(dims[1,])
          else
            mats<-1:NCOL(dims)
          dims<-dims[2,]         #Get ncols
        }
        if((!any(is.null(dims)))&&(length(unique(dims))==1)&&(all(mats>0))){
          out<-array(dim=c(mats[length(mats)],dims[1],dims[1]))
          for(i in 1:length(mats))
            out[(c(0,mats)[i]+1):(mats[i]),,]<-g[[i]]
        }else
          out<-g
      }else
        out<-g
    }
  }else{
    if((!simplify)&&(length(dim(g))==3)){  #Coerce to a list of matrices?
      out<-list()
      for(i in 1:dim(g)[1])
        out[[i]]<-g[i,,]
    }else
      out<-g
  }
  #Return the result
  out
}
