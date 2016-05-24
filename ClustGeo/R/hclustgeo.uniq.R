#' @export
#' @name hclustgeo.uniq
#' @title Hierarchical clustering based on geographical distance and on distance between individuals mesured on numerical variables.
#' @description This version optimize a global criterion of quality of classes H= alpha I_intra_var + (1-alpha) I_intra_geo
#' @param data a data frame with \code{n} rows and \code{p} colums containing measures of individuals on variables
#' @param D.geo a matrix of size \code{n} by \code{n} containing geographical distances between individuals
#' @param alpha a real between 0 and 1 to fit the global criterion to optimize
#' @return {clust}{an object of class hclustgeo.uniq containing the result of hclustgeo}
#' @keywords internal

hclustgeo.uniq <- function(data, D.geo, alpha) {
  
  cl <- match.call()
  
  n<-nrow(data)
  name <- rownames(data)
  
  
  if (nrow(D.geo)!=ncol(D.geo))
    stop("\"D.geo\" must be a symmetric matrix of dimension n")
  
  if(length(which((rownames(data)==rownames(D.geo))==FALSE))!=0)
    stop("rownames of \"data\" must be in the same order as rownames of \"D.geo\"")

  
  
  MAXVAL <- 1.0e12
  
  
  #Weights of individuals all equals to 1/n
  wt<-rep(1/n,n)
  wt.init<-wt

  #D is the matrix of euclid distances calculated on 'data'
  D<-as.matrix(dist(data,upper=TRUE,diag=TRUE))
  D.red<-D/max(D)
  D.red.2<-D.red^2
  #diss.data is the scaled matrix of aggreg measure based on D (so based on data)
  diss.data<-(1/(2*n))*D.red.2
  #Save of the no scaled distance matrix (no squared) on data
  dist.data.save<-D
  
  #D.geo 
  dist.geo.save<-D.geo  
  D.geo.red<-D.geo/max(D.geo)
  D.geo.red.2<-D.geo.red^2
  #diss.data is the scaled matrix of aggreg measure based on D.geo 
  diss.geo<-(1/(2*n))*D.geo.red.2
  
  #diss is the "total" matrix of aggreg measures
  diss<-alpha*diss.data + (1-alpha)*diss.geo
  
  
  flag <- rep(1, n)                          # active/dead indicator
  a <- rep(0, n-1)                           # left subnode on clustering
  b <- rep(0, n-1)                           # right subnode on clustering
  ia <- rep(0, n-1)                          # R-compatible version of a
  ib <- rep(0, n-1)                          # R-compatible version of b
  lev <- rep(0, n-1)                         # level or criterion values
  card <- rep(1, n)                          # cardinalities
  mass <- wt
  order <- rep(0, n)                         # R-compatible order for plotting
  
  
  
  
  nnsnnsdiss <- getnns(diss, flag)           # call to function getnns
  clusmat <- matrix(0, n, n)                 # cluster memberships
  for (i in 1:n) clusmat[i,n] <- i           # init. trivial partition
  
  
  #ncl=2
  seq.print<-seq(from=1,to=n-1,by=50)
  for (ncl in (n-1):1) {                      # main loop 
    # check for agglomerable pair
    if(abs(ncl-n)%in%seq.print)
      cat ("iter=",abs(ncl-n+1)," /",n-1,"\n",sep="")
    #cat("ncl=",ncl,"\n")
    #cat("clusmat=","\n")
    #print(clusmat)
    #cat("dissim=","\n")
    #print(diss)
    minobs <- -1;  
    mindis <- MAXVAL;
    for (i in 1:n) {
      if (flag[i] == 1) {
        if (nnsnnsdiss$nndiss[i] < mindis) {
          mindis <- nnsnnsdiss$nndiss[i]
          minobs <- i
        }
      }
    }
    #cat("mindis=",mindis,"\n")
    #cat("minobs=",minobs,"\n")
    
    # find agglomerands clus1 and clus2, with former < latter
    if (minobs < nnsnnsdiss$nn[minobs]) {
      clus1 <- minobs
      clus2 <- nnsnnsdiss$nn[minobs]
    }
    if (minobs > nnsnnsdiss$nn[minobs]) {
      clus2 <- minobs
      clus1 <- nnsnnsdiss$nn[minobs]
    }
    #cat("clus1=",clus1,"\n")
    #cat("clus2=",clus2,"\n")
    
    # So, agglomeration of pair clus1 < clus2 defines cluster ncl
    
    #------------------------------------ Block for subnode labels 
    a[ncl] <- clus1                       # aine, or left child node
    b[ncl] <- clus2                       # benjamin, or right child node
    # Now build up ia, ib as version of a, b which is R-compliant
    
    #cat("card[clus1] ==",card[clus1],"\n")
    #cat("card[clus2] ==",card[clus2],"\n")
    
    if (card[clus1] == 1) ia[ncl] <- (-clus1)     # singleton
    if (card[clus2] == 1) ib[ncl] <- (-clus2)     # singleton
    if (card[clus1] > 1) {                # left child is non-singleton
      lastind <- 0
      for (i2 in (n-1):(ncl+1)) {        # Must have n-1 >= ncl+1 here
        if (a[i2] == clus1) lastind <- i2    # Only concerns a[i2]
      }
      ia[ncl] <- n - lastind             # label of non-singleton
    }
    if (card[clus2] > 1) {                # right child is non-singleton
      lastind <- 0
      for (i2 in (n-1):(ncl+1)) {        # Must have n-1 >= ncl+1 here
        if (a[i2] == clus2) lastind <- i2    # Can only concern a[i2]
      }
      ib[ncl] <- n - lastind             # label of non-singleton
    }
    if (ia[ncl] > 0 || ib[ncl] > 0) {     # Check that left < right
      left <- min(ia[ncl],ib[ncl])
      right <- max(ia[ncl],ib[ncl])
      ia[ncl] <- left                    # Just get left < right
      ib[ncl] <- right
    }
    #--------------------------------------------------------------------
    
    lev[ncl] <- mindis
    for (i in 1:n) {
      clusmat[i,ncl] <- clusmat[i,ncl+1]
      if (clusmat[i,ncl] == clus2) clusmat[i,ncl] <- clus1
    }
    # Next we need to update diss.geo and diss.data array
    for (i in 1:n) {
      if ( (i != clus1) && (i != clus2) && (flag[i] == 1) ) {
        
        diss.data[clus1,i] <- 
          ((mass[clus1]+mass[i])/(mass[clus1]+mass[clus2]+mass[i]))*diss.data[clus1,i] +
          ((mass[clus2]+mass[i])/(mass[clus1]+mass[clus2]+mass[i]))*diss.data[clus2,i] -
          (mass[i]/(mass[clus1]+mass[clus2]+mass[i]))*diss.data[clus1,clus2] 
        diss.data[i,clus1] <- diss.data[clus1,i]
        
        
        diss.geo[clus1,i] <- 
          ((mass[clus1]+mass[i])/(mass[clus1]+mass[clus2]+mass[i]))*diss.geo[clus1,i] +
          ((mass[clus2]+mass[i])/(mass[clus1]+mass[clus2]+mass[i]))*diss.geo[clus2,i] -
          (mass[i]/(mass[clus1]+mass[clus2]+mass[i]))*diss.geo[clus1,clus2] 
        diss.geo[i,clus1] <- diss.geo[clus1,i]
        
        
      }
    }
    mass[clus1] <- mass[clus1] + mass[clus2]    # Update mass of new cluster
    card[clus1] <- card[clus1] + card[clus2]    # Update card of new cluster
    # Cluster label clus2 is knocked out; following not nec. but no harm
    flag[clus2] <- 0
    nnsnnsdiss$nndiss[clus2] <- MAXVAL
    mass[clus2] <- 0.0
    
    for (i in 1:n) {
      diss.data[clus2,i] <- MAXVAL
      diss.data[i,clus2] <- diss.data[clus2,i]
      
      diss.geo[clus2,i] <- MAXVAL
      diss.geo[i,clus2] <- diss.geo[clus2,i]
    }
    
    # Finally update nnsnnsdiss$nn and nnsnnsdiss$nndiss
    # i.e. nearest neighbors and the nearest neigh. dissimilarity
    diss<-alpha*diss.data + (1-alpha)*diss.geo
    
    nnsnnsdiss <- getnns(diss, flag)
    
  }
  cat ("iter=",n-1," /",n-1,"\n",sep="")
  
  temp <- cbind(a,b)
  merge2 <- temp[nrow(temp):1, ]
  temp <- cbind(ia,ib)
  merge <- temp[nrow(temp):1,]
  dimnames(merge) <- NULL
  # merge is R-compliant; later suppress merge2
  
  #-------------------------------- Build R-compatible order from ia, ib
  orderlist <- c(merge[n-1,1], merge[n-1,2])
  norderlist <- 2
  for (i in 1:(n-2)) {           # For precisely n-2 further node expansions
    #cat("i=",i,"\n")
    for (i2 in 1:norderlist) {       # Scan orderlist
      #cat("i2=",i2,"\n")
      if (orderlist[i2] > 0) {     # Non-singleton to be expanded
        tobeexp <- orderlist[i2]
        if (i2 == 1) {
          orderlist <- c(merge[tobeexp,1],merge[tobeexp,2],
                         orderlist[2:norderlist])
        }
        if (i2 == norderlist) {
          orderlist <- c(orderlist[1:(norderlist-1)],
                         merge[tobeexp,1],merge[tobeexp,2])
        }
        if (i2 > 1 && i2 < norderlist) {
          orderlist <- c(orderlist[1:(i2-1)], 
                         merge[tobeexp,1],merge[tobeexp,2],
                         orderlist[(i2+1):norderlist])
        }
        norderlist <- length(orderlist)
      }
    }
  }
  orderlist <- (-orderlist)
  class(orderlist) <- "integer"
  
  xcall <- "hierclust(a,wt)"
  class(xcall) <- "call"
  #clusmat=clusmat
  #labels=as.character(1:n)
  
  retlist <- list(merge=merge,height=as.single(lev[(n-1):1]),order=orderlist,
                  labels=dimnames(a)[[1]],method="minvar",call=xcall,
                  dist.method="euclidean-factor")
  retlist <- list(call=cl, merge=merge,height=lev[(n-1):1],order=orderlist,alpha=alpha, dist.var=dist.data.save, dist.geo=dist.geo.save, 
                  wt=wt.init, data=data)
  class(retlist) <- "hclustgeo.uniq"
  retlist
}