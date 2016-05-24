netdist <- function(x, ...) UseMethod("netdist")

## Method netdist for matrix
netdist.matrix <- function(x, h=NULL, d="HIM", ga=NULL, components=TRUE, ...){

  if (is.null(h))
    stop("Need to provide a second matrix to compute the distance")
  
  DISTANCE <- c("HIM","IM","H",
                "ipsen","Ipsen","IpsenMikhailov","Ipsen-Mikhailov",
                "hamming","Hamming")
  d <- pmatch(d, DISTANCE)
  if(d==2L | (d>=4  & d<=7L))
    d <- 2L
  if(d==3L |(d>=8L & d<=9L))
    d <- 3L
  
  if(is.na(d))
    stop("invalid distance", call. =FALSE)
  if(d == -1)
    stop("ambiguous distance", call. =FALSE)
  
  Call <- match.call()
  
  ##add a check so that an unexisting parameter cannot be passed
  id.Call <- match( names(Call),c("x", "h", "d", "ga","components","n.cores","verbose", "rho"), nomatch=0)
  if(sum(id.Call[-1]==0)==1){
    warning("For netdist function the parameter '",names(Call)[which(id.Call==0)[2]],"' will be ignored",call.=FALSE)
  }
  if(sum(id.Call[-1]==0)>1){
    msg <- "For netdist function, the following parameters will be ignored:\n"
    for(i in which(id.Call==0)[-1]){
      msg <- paste(msg,"'",names(Call)[i],"'\n",sep="")
    }
    warning(msg,call.=FALSE)
  }
  if(is.na(d))
    stop("invalid distance")
  if(d == -1)
    stop("ambiguous distance")
  
  ##check on components argument
  if(is.null(components)){
    if(d==1){
      comp <- TRUE
    }else{
      comp <- FALSE
    }
  }else{
    comp <- components
    if(d==1){
      if(!is.logical(comp))
        stop("components must be TRUE or FALSE")
    }else{
      comp <- FALSE
      # warning("components parameter will be ignored", call. = FALSE)
    }
  }

  ##check on ga passing through ipsen function
  if(is.null(ga)){
    if(d==2){
      warning("The ga parameter will be automatically defined.", call.=FALSE)
    }
  }else{
    if(!is.numeric(ga) && !is.null(ga))
      stop("ga must be numeric",call.=FALSE)
  }
  
  ## Create two classes of the same type.
  ## One for each input
  g1 <- g2adj(x)
  g2 <- g2adj(h)
  
  ## Check for dimension and directionality of g and h
  ## Return two different error messages
  if ((g1$N == g2$N)){
    if((g1$tag == g2$tag)){
      ## Create the list of adjacency matrix
      myadj <- list(d=DISTANCE[d],G=list(g1$adj,g2$adj),N=g1$N,tag=g1$tag)
    } else {
      stop("Not conformable graph g and h: one is directed while the other is undirected", call.=FALSE)
    }
  } else {
    stop("Not conformable graph g and h: they have different dimensions", call.=FALSE)
  }
  
  ## Check the distance chosen and compute it
  ## HIM
  if (myadj$d=="HIM"){
    mylap <- list(L=list(Lap(g1$adj), Lap(g2$adj)),N=g1$N, tag=g1$tag)
    dd <- him(list(ADJ=myadj,LAP=mylap),  ga=ga,  components=comp, ltag=FALSE, ...)
  }
  ## IM
  if (myadj$d=="IM"){
    mylap <- list(L=list(Lap(g1$adj), Lap(g2$adj)),N=g1$N, tag=g1$tag)
    ## mylap <- list(L1=Lap(g1$adj),L2=Lap(g2$adj),N=g1$N, tag=g1$tag)
    dd <- ipsen(mylap,ga=ga, ...)
  }
  ## H
  if (myadj$d=="H"){
    dd <- hamming(myadj)
  }

  ## Return distance list
  return(dd)
}

## Register S3 method for object Matrix
netdist.Matrix <- netdist.matrix

## Register S4 methods for matrix, Matrix and data.frame objects
setMethod("netdist","matrix", netdist.matrix)
setMethod("netdist","Matrix", netdist.matrix)
setMethod("netdist","data.frame", netdist.matrix)


## Register S3 and S4 methods for igraph objects is igraph is available
netdist.igraph <- netdist.matrix
setMethod("netdist","igraph",netdist.matrix)


## Method netdist for list of adjacency matrices
##--------------------------------------------------
## Implementation for kernel distance
## Output a distance matrix
netdist.list <- function(x, d="HIM", ga=NULL, components=TRUE, ...){
  DISTANCE <- c("HIM","IM","H",
                "ipsen","Ipsen","IpsenMikhailov","Ipsen-Mikhailov",
                "hamming","Hamming")
  d <- pmatch(d, DISTANCE)
  if(d==2L | (d>=4  & d<=7L))
    d <- 2L
  if(d==3L |(d>=8L & d<=9L))
    d <- 3L

  ## Check distance type
  if(is.na(d))
    stop("invalid distance", call. =FALSE)
  if(d == -1)
    stop("ambiguous distance", call. =FALSE)

  ## Store the function call
  Call <- match.call()
  id.Call <- match( names(Call),c("x", "d", "ga","components", "n.cores", "verbose", "rho"), nomatch=0)

  ## add a check so that an unexisting parameter cannot be passed
  if(sum(id.Call[-1]==0)==1){
    warning("The parameter '",names(Call)[which(id.Call==0)[2]],"' will be ignored",call.=FALSE)
  }
  if(sum(id.Call[-1]==0)>1){
    msg <- "The following parameters will be ignored:\n"
    for(i in which(id.Call==0)[-1]){
      msg <- paste(msg,"'",names(Call)[i],"'\n",sep="")
    }
    warning(msg,call.=FALSE)
  }
  
  ##check if need to return all components
  if(is.null(components)){
    if(d==1){
      comp <- TRUE
    } else {
      comp <- FALSE
    }
  } else {
    comp <- components
    if(d==1){
      if(!is.logical(comp))
        stop("components must be TRUE or FALSE")
    } else {
      comp <- FALSE
      warning("components parameter will be ignored", call. = FALSE)
    }
  }
  
  ##check on ga pass through ipsen function
  if(is.null(ga)){
    if(d==2){
      warning("The ga parameter will be automatically defined.", call.=FALSE)
    }
  }else{
    if(!is.numeric(ga) && !is.null(ga))
      stop("ga must be numeric",call.=FALSE)
  }


  ## Distance computation
  ##------------------------------
  ## Create the class
  tmp <- lapply(x, g2adj)
  if (d == 1L | d == 2L)
    laplist <- list()
  
  adjlist <- list()
  N <- tmp[[1]]$N
  tag <- tmp[[1]]$tag

  ## Create a list of adjacency matrices
  for (i in 1:length(tmp)){
    g <- tmp[[i]]
    if (g$N == N && g$tag == tag){
        if (d == 1L | d == 2L)
          laplist[[i]] <- Lap(g$adj)
        adjlist[[i]] <- g$adj
    } else {
      stop(paste("Element",i,"not of the same length!"), call.=FALSE)
    }
    myadj <- list(d=DISTANCE[d],G=adjlist,N=N,tag=tag)
  }
  
  ## Compute the distance

  ## HIM
  if (myadj$d=="HIM"){
    mylap <- list(L=laplist,N=N, tag=tag)
    dd <- him(list(ADJ=myadj,LAP=mylap), ga=ga,  components=comp, ltag=TRUE, ...)
  }

  ## IM
  if (myadj$d=="IM"){
    mylap <- list(L=laplist,N=N, tag=tag)
    dd <- ipsen(mylap,ga=ga, ...)
    dd <- list(IM=dd)
  }

  ## H
  if (myadj$d=="H"){
    dd <- hamming(myadj)
    dd <- list(H=dd)
  }
  
  ## Give names to the matrices
  if (!is.null(names(x))){
    ## Create a list from the matrix
    if (is.list(dd)){
      dd <- lapply(dd,function(y,x){colnames(y) <- rownames(y) <- names(x)
                                    return(y)}, x=x)
    } else {
      colnames(dd) <- rownames(dd) <- names(x)
      dd <- list(HIM=dd)
    }
  } else {
    ## Return names in the list if HIM with no components has been set
    if (is.matrix(dd) && d==1L){
      dd <- list(HIM=dd)
    } 
  }

  ## Return a distance matrix
  return(dd)
}

## Register S4 method for lists
setMethod("netdist", "list", netdist.list)


## Prepare the matrix for computing distance if the graph is directed
transfmat <- function(x){
  ## Check if the x matrix is symmetric (undirected graph)
  ## Otherwise it returns a list with a matrix like:
  ## |zeros    t(A)|
  ## |  A     zeros|
  ##
  Adj <- x
  n <- ncol(Adj)
  tag <- "undir"
  
  ## If the graph is directed create a new matrix 
  if (!isSymmetric(x,check.attributes=FALSE, check.names=FALSE)){
    zero <- matrix(0, nrow=n, ncol=n)
    tmp <- Matrix::rBind(Matrix::cBind(zero,t(Adj)),Matrix::cBind(Adj,zero))
    Adj <- tmp
    tag <- "dir"
  }

  ## Normalize the weights
  if (any(Adj > 1) || any(Adj < 0)){
    warning("Edge weight should be >= 0 and <= 1, scaling has been automatically applied!", call.=FALSE)
    Adj <- (Adj - min(Adj)) / (max(Adj) - min(Adj))
  }
  return(list(adj = Adj, tag = tag, N=n))
}

## Create the adjacency matrix structure
##----------------------------------------
g2adj <- function(x,...) UseMethod("g2adj")

g2adj.igraph <- function(x,...,type="both"){
  if (!is.null(get.edge.attribute(x,"weight"))){
    WW <- "weight"
  } else {
    warning("No weight attribute to the graph object\nCompute binary adjacency matrix", call.=FALSE)
    WW <- NULL
  }
  Adj <- get.adjacency(x,type=type,attr=WW,sparse=TRUE)
  diag(Adj) <- 0
  ll <- transfmat(Adj)

  return(ll)
}
setMethod("g2adj","igraph",g2adj.igraph)

g2adj.matrix <- function(x,...){
  ll <- transfmat(x)
  return(ll)
}
setMethod("g2adj","matrix",g2adj.matrix)
setMethod("g2adj","Matrix",g2adj.matrix)

g2adj.data.frame <- function(x, ...){
  x <- apply(x,2,as.numeric)
  ll <- transfmat(x)
  return(ll)
}
setMethod("g2adj","data.frame",g2adj.data.frame)

## Generical Laplacian
##----------------------------------------
Lap <- function(x,...){
  D <- apply(x,2,sum)
  L <- -x
  diag(L) <- D
  return(L)
}

## Him distance
##----------------------------------------
him <- function(object,ga=NULL, components=TRUE, ltag=FALSE, rho=1, ...){
  ipd <- ipsen(object$LAP, ga, ...)
  had <- hamming(object$ADJ)
  gloc <- sqrt(1/(1+rho)) * sqrt(had**2+ rho*(ipd**2))
  if (length(gloc)==1)
    names(gloc) <- "HIM"
  if(components==TRUE){
    if (ltag){
      dist <- list(H=had,IM=ipd,HIM=gloc)
    } else {
      dist <- c(had, ipd, gloc)
      names(dist) <- c("H","IM","HIM")
    }
  } else {
    dist <- gloc
    if (!ltag)
      names(dist) <- "HIM"
  }
  return(dist)
}
