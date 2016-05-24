#Given an nxp matrix m and a function f,
# returns the pxp matrix got by applying f  to all pairs of columns of  m.

colpairs <- function(m,f,diag=0,na.omit=FALSE,...){
    flocal <- function(i,j) 
       if (!is.null(diag) && (i == j))
           diag
        else {
	   x <- m[,i]
	   y <- m[,j]
	   if (na.omit) {
              d <- na.omit(cbind(x,y))
	      x <- d[,1]
	      y <- d[,2]}
	   f(x,y,...) }
       
   p <- ncol(m)
   m1 <- matrix(rep(1:p,p),nrow=p,ncol=p)
   ind <- mapply("c",m1,t(m1))
   ans <- apply(ind,2, function(i) flocal(i[1],i[2]))
   ans <- matrix(ans,nrow=p,ncol=p)
   colnames(ans) <- colnames(m)
   rownames(ans) <- colnames(m)
   ans
   }
	

km2 <- function(x,y){
   x <- x - mean(x)
   y <- y - mean(y)
   sum(x*x)+ sum(y*y)
   }

# Computes the sum of all distances between pairs of
# objects whose coordinates are contained in x and y.
gtot <- function(x,y,...)
	2*sum(dist(cbind(x,y),...))


# Computes the average total  distance from one object to all other 
# objects, where x and y contain the object cordinates.

gave <- function(x,y,...)
   2*sum(dist(cbind(x,y),...))/length(x)


# Computes the cluster diameter- the maximum distance between
# objects whose coordinates are contained in x and y.

diameter <- function(x,y,...){
   d <- dist(cbind(x,y),...)
   max(d)
}

# Computes the cluster star distance- the minimum of the total distance from
# one object to another, where x and y contain the object cordinates.

star <- function(x,y,...){
   d <- vec2distm(dist(cbind(x,y),...))
   min(apply(d,2,sum))
   }


# Computes the silhouette distance of a partition of the objects in
# x and y, where group contains the object memberships.

sil <- function(x,y,groups,...){
  # require(cluster)
   igroups <- unclass(factor(groups))
   d <- dist(cbind(x,y),...)
   s <- silhouette(igroups,d)
   summary(s)$avg.width
}

# Computes the agglomerative coefficient, from agnes.

ac <- function(x,y,...){
#   require(cluster)
   ag <- agnes(cbind(x,y),keep.diss=FALSE,keep.data=FALSE,...)
   ag$ac
}

# Computes the total line length in a parallel coordinate plot
# of x and y.
pclen <- function(x,y) sum(abs(y-x))

# Computes the average (per object) line length in a parallel coordinate plot
# where each x object is connected to all y objects.
pcglen <- function(x,y)
   sum(outer(x,y,function(a,b) abs(a-b)))/length(x)
	

# Applies the function gfun  to each group of x and y values
# and combines the results using the function cfun.
#(...) arguments are passed to gfun.

partition.crit <- function(x,y,groups,gfun= gave,cfun=sum,...){
   dgroups <- unique(groups)
   gm <- sapply(dgroups,function(g) gfun(x[groups==g],y[groups==g],...))
   cfun(gm)
}





