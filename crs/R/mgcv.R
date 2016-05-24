
#############################################################
## Type 1 tensor product methods start here (i.e. Wood, 2006)
#############################################################

tensor.prod.model.matrix<-function(X)
# X is a list of model matrices, from which a tensor product model matrix is to be produced.
# e.g. ith row is basically X[[1]][i,]%x%X[[2]][i,]%x%X[[3]][i,], but this routine works 
# column-wise, for efficiency
{ m<-length(X)
  X1<-X[[m]]
  n<-nrow(X1)
  if (m>1) for (i in (m-1):1)
  { X0<-X1;X1<-matrix(0,n,0)
    for (j in 1:ncol(X[[i]]))
    X1<-cbind(X1,X[[i]][,j]*X0)
  }
  X1
} ## end tensor.prod.model.matrix

uniquecombs<-function(x) {
## takes matrix x and counts up unique rows
## `unique' now does this in R
if (is.null(x)) stop("x is null")
if (is.null(nrow(x))) stop("x has no row attribute")
if (is.null(ncol(x))) stop("x has no col attribute")
ind <- rep(0,nrow(x))
res<-.C("RuniqueCombs",x=as.double(x),ind=as.integer(ind),
        r=as.integer(nrow(x)),c=as.integer(ncol(x)))
n <- res$r*res$c
x <- matrix(res$x[1:n],res$r,res$c)
attr(x,"index") <- res$ind+1 ## C to R index gotcha
x
}

