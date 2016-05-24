
iplot <- function(x,...){
  UseMethod("iplot")
}

iplot.graphNEL <- function(x,...){
  ig <- igraph::igraph.from.graphNEL(x)
  igraph::V(ig)$label <- igraph::V(ig)$name
  igraph::V(ig)$size  <- 50
  ig$cex   <-  4
                                        #ig$layout   <- layout.graphopt
                                        #ig$layout <- layout.kamada.kawai
  ig$layout <- igraph::layout.lgl
  plot(ig,
       vertex.label.family="Helvetica",
       edge.label.family="Helvetica",
       vertex.label.cex=2,
       edge.label.cex=2)
}



.dgCMatrix <- function(data=NA, nrow=1, ncol=1, byrow=FALSE, dimnames=NULL,
                      sparse = TRUE, doDiag = TRUE, forceCheck = FALSE){
  as(Matrix(data=data, nrow=nrow, ncol=ncol, dimnames=dimnames, sparse=TRUE), "dgCMatrix")
}


### rowmat2list and colmat2list:
### ----------------------------
## Turns a matrix into a list, either by row or by column.
## Notice: finding unique rows in a matrix can be speeded up this way.

matrix2list <- function(x, byrow=TRUE){
  if (byrow)
    rowmat2list(x) # cpp implementation
  else
    colmat2list(x) # cpp implementation
}

## Returns matrix n x 2 matrix with indices of non-zero
## entries in matrix m
## FIXME: which.arr.ind: Fails on sparse matrices!!
## FIXME: -> remove after check downstram!!
## FIXME: -> which_matrix_index is Cpp implementation
which.arr.ind<-function(m){
  nr  <- nrow(m)
  nc  <- ncol(m)
  rr <- rep.int(1:nr, nc)
  cc <- rep(1:nc, each=nr)
  cbind(rr[m!=0L], cc[m!=0L])
}


## lapplyMatch: same as but much faster than
## lapply(xlist, function(gg) match(gg, set))
##
lapplyV2I <- lapplyMatch <- function(xlist, set){lapply(xlist, function(gg) match(gg, set))}

## lapplyI2C: same as but faster than
## lapply(xlist, function(x) set[x])
lapplyI2V <- function (xlist, set) {lapply(xlist, function(xx) set[xx])}

##
## Calculate logL for N(0,\Sigma) model.
##
## Sigma = Covariance matrix parameter
## K     = Sigma inverse
## S     = sample covariance matrix
## n     = sample size
##
ell <- function(Sigma, S, n){

  shdet <- function(Sigma){
    prod(eigen(Sigma)[[1]])
  }
  p <- dim(S)[1]
  const <- -n*p/2*log(2*pi)
  return(const-n/2*log(shdet(Sigma))
         -n/2*sum(diag( solve(Sigma)%*%S )) )
}

ellK <- function (K, S, n)
{
    value <- (n/2) * (log(det(K)) - sum(rowSums(K * S)))
    return(value)
}

cov2pcor <- function(V){
  ans <- -cov2cor(solve(V))
  diag(ans) <- -diag(ans)
  ans
  }

conc2pcor <- function(K){
  ans <- -cov2cor(K)
  diag(ans)<-1
  ans
}




## Codes a p x 2 matrix of characters or a list with pairs
## of characters into a vector of numbers.

## FIXME: pairs2num: Cpp implementation
pairs2num <- function(x, vn, sort=TRUE){
  if (class(x)!="matrix"){
    if (is.null(x))
      return(NULL)

    if (inherits(x,"list"))
      x <- do.call(rbind,x)
    else {
      if (inherits(x,"character"))
        x <- matrix(x,nrow=1)
    }
  }
  # From here x should be a p x 2 matrix

  dd <- dim(x)
  if (dd[1L]==0){
      return(numeric(0))
  } else {
      if (sort){
          i     <- x[,2L]< x[,1L]
          c1    <- i+1L
          c2    <- -1L*(i-1L) + 1L
          x  <- cbind(
                      x[cbind(seq_along(c1),c1)],
                      x[cbind(seq_along(c2),c2)])
        }
      ans       <- match(x,vn)
      dim(ans)  <- dim(x)
      colSumsPrim(t.default(ans) * c(100000,1))
      ## ans[,1L] <- ans[,1L] * 100000L
##       rowSumsPrim(ans)
    }
}

####################################################
####
#### Create all possible pairs from a vector
#### or all possible pairs combining one element
#### from each of two vectors
####
#### NOTICE: If y is not NULL then x and y must be disjoint
####
####################################################

## FIXME: names2pairs: Cpp implementation (code is in mail somewhere)
names2pairs <- function(x, y=NULL, sort=TRUE, result="list"){
  result <- match.arg(result, c("list","matrix"))
  lenx <- length(x)
  leny <- length(y)
  if (leny==0){
    if (lenx==1){
      if (result=="matrix")
        return(matrix(nrow=0,ncol=2))
      else
        return(list())
    } else {
      cc   <- combnPrim(1:length(x),2)
      ans  <- x[cc]
      dim(ans) <- dim(cc)
      if (sort){
        idx <- ans[1,]>ans[2,]
        ans[1:2,idx] <- ans[2:1,idx]
      }
      if (result=="matrix")
        return(t.default(ans))
      else
        return(colmat2list(ans))
    }
  } else {
    ans <- cbind(rep(x, each=leny), rep(y, times=lenx))
    if (sort){
      idx <- ans[,1]>ans[,2]
      ans[idx, 1:2] <- ans[idx,2:1]
    }
    if (result=="matrix")
      return(ans)
    else
      rowmat2list(ans)
  }
}

