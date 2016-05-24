###Now for %*%
setClass("SparseplusLowRank",representation(x="sparseMatrix",a="matrix",b="matrix"))
splr <- function(x,a=NULL,b=NULL){
  # x+ab'
  dx=dim(x)
  if(is.null(a))b=NULL
  if(is.null(b))a=NULL
  if(!is.null(a)){   
    da=dim(a)
    db=dim(b)
    if(da[1]!=dx[1])stop("number of rows of x not equal to number of rows of a")
    if(db[1]!=dx[2])stop("number of columns of x not equal to number of rows of b")
    if(da[2]!=db[2])stop("number of columns of a not equal to number of columns of b")
  }
  new("SparseplusLowRank",x=x,a=a,b=b)
}

.leftmult=function(x,y){
  #y is splr, x is matrix
  a=y@a
  b=y@b
  sx=y@x
  if(is.null(a)|is.null(b))x %*% sx
  else{
    part1=x %*% sx
    part2=x%*% a
    part2=part2 %*% t(b)
    part1+part2
  }
}
.rightmult=function(x,y){
  #x is splr, y is matrix
  a=x@a
  b=x@b
  sx=x@x
  if(is.null(a)|is.null(b))sx%*%y
  else{
    part1=sx %*% y
    part2=t(b) %*% y
    part2=a%*% part2
    part1+part2
  }
}
Frobsmlr=function(x,a,b,nx=NULL){
###Note: this computes the Frobenius norm of the Difference between x and ab'
if(is.null(nx))nx=norm(x,type="F")

  xab=as.matrix(x%*%b)
  xab=sum(xab*a)
  aa=t(a)%*%a
  bb=t(b)%*%b
  ab=sum(aa*bb)
  sqrt(pmax(0,nx^2-2*xab+ab))
}

Frob=function(Uold,Dsqold,Vold,U,Dsq,V){
denom=sum(Dsqold^2)
utu=Dsq* (t(U)%*%Uold)
vtv=Dsqold* (t(Vold)%*%V)
uvprod= sum(diag(utu%*%vtv))
num=denom+sum(Dsq^2) -2*uvprod
num/max(denom,1e-9)
}
  

  


setMethod("%*%",signature(x="SparseplusLowRank",y="Matrix"),.rightmult)
setMethod("%*%",signature(x="Matrix",y="SparseplusLowRank"),.leftmult)
setMethod("%*%",signature(x="SparseplusLowRank",y="ANY"),.rightmult)
setMethod("%*%",signature(x="ANY",y="SparseplusLowRank"),.leftmult)
setMethod("dim", signature(x = "SparseplusLowRank"),
	  function(x) dim(x@x), valueClass = "integer")
setMethod("norm",signature(x="SparseplusLowRank",type="character"),
          function(x,type,...){
            switch(type,
                   "F"=Frobsmlr(x=x@x,a=x@a,b=-x@b),
                   stop("invalid 'type'"))
          },valueClass="numeric")


### rowSums, colSums, rowMeans, colMeans 
.rsum=function(x,...){
  #x is SparseplusLowRank matrix
  rx=rowSums(x@x)
  cb=colSums(x@b)
  drop(rx+x@a%*%cb)
}
.csum=function(x,...){
  #x is SparseplusLowRank matrix
  cx=colSums(x@x)
  ca=colSums(x@a)
 drop( cx+x@b%*%ca)
}
.rmean=function(x,...){
  #x is SparseplusLowRank matrix
  rx=rowMeans(x@x)
  cb=colMeans(x@b)
  drop(rx+x@a%*%cb)
}
.cmean=function(x,...){
  #x is SparseplusLowRank matrix
  cx=colMeans(x@x)
  ca=colMeans(x@a)
  drop(cx+x@b%*%ca)
}
as.matrix.splr=function(x,...)  as.matrix(x@x)+x@a%*%t(x@b)
setMethod("rowSums","SparseplusLowRank",.rsum)
setMethod("colSums","SparseplusLowRank",.csum)
setMethod("rowMeans","SparseplusLowRank",.rmean)
setMethod("colMeans","SparseplusLowRank",.cmean)
setMethod("as.matrix","SparseplusLowRank",as.matrix.splr)

colsum2.splr <-
  function(z){
###Compute column sums of (X+AB^T)^2  (element wise squares)
### where z is stored as SparseplusLowRank
    X=z@x
    A=z@a
    B=z@b
    X2=X
    X2@x=X@x^2
    s1=colSums(X2)
    s2=colSums( as.matrix(t(A)%*%X)*t(B))

    c2=t(A)%*%A
    c2=B%*%c2
    s3=rowSums(B*c2)
    s1+2*s2+s3
  }


rowsum2.splr <-
  function(z){
###Compute row sums of (X+AB^T)^2  (element wise squares)
### where z is stored as SparseplusLowRank

    X=z@x
    A=z@a
    B=z@b
    X2=X
    X2@x=X@x^2
    s1=rowSums(X2)
    s2=rowSums( as.matrix(X%*%B)*A)

    c2=t(B)%*%B
    c2=A%*%c2
    s3=rowSums(A*c2)
    s1+2*s2+s3
  }

