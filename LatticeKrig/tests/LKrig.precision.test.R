
# LatticeKrig
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library( LatticeKrig)
options( echo=FALSE)

##########################################
  test.for.zero.flag<- 1

# test building MRF
  mx<- rbind(c(5,6))
 m<- mx[1,1]* mx[1,2]
temp<- matrix( 1:m, mx[1,1],mx[1,2])
 m1<- mx[1,1]
 m2<- mx[1,2]
   I.c<- temp 
   I.B<-  cbind(      rep(NA,m1),           temp[,1:(m2-1)])
   I.T<-  cbind(       temp[,2:m2],        rep(NA,m1))
   I.L<-  rbind(     rep(NA , m2), temp[ 1:(m1-1),] )
   I.R<-  rbind(  temp[2:m1 , ],          rep( NA, m2))

   Bi<- rep( 1:30, 5)
   Bj<- cbind( c(I.c), c(I.T), c(I.B), c( I.L), c(I.R))
   atest<- matrix( (1:m)*5, m1,m2)
   values<- cbind( c(atest), rep(-1,m),rep(-1,m),rep(-1,m),rep(-1,m) )
   good<- !is.na(c(Bj))
   Bi<- Bi[good]
   Bj<- Bj[good]
   values<- c(values)[c(good)]
   LKinfo<- LKrigSetup(  cbind( c( 1,5), c( 1,6)), NC=6, a.wght=matrix(atest, 5,6), nlevel=1, alpha=1, NC.buffer=0)
   obj<- LKrigSAR( LKinfo, Level=1)

   test.for.zero( cbind( Bi,Bj), obj$ind, tag="MRF index")
   test.for.zero( values, obj$ra, tag="MRF value")

# 

   LKinfo<- LKrigSetup(  cbind( c( 1,5), c( 1,6)), NC=6, a.wght=matrix(atest, 5,6),alpha=1, nlevel=1, NC.buffer=0)
   obj<- LKrigSAR( LKinfo, Level=1)
   obj2<-spind2full( obj)
   test2<- matrix( obj2[8,], 5,6)
   test0<- matrix( 0, 5,6)
   test0[3,2]<- 5*8
   test0[2,2]<- test0[3,1]<- test0[4,2]<- test0[3,3] <- -1
   test.for.zero( test2, test0, tag="MRF weight placement")


############ end of simple MRF tests


# now everything
# TEST INDEXING
#
  set.seed(123)
  x1<- cbind( runif( 10), runif(10))
  LKinfo<- LKrigSetup( cbind( c( -1,1), c( -1,1) ),
             nlevel=3, NC=4, a.wght=c(5,6,7), alpha=c(6,6,6) )
  X<- LKrig.basis(x1, LKinfo)
  X<- as.matrix(X)
# check on normalization to unit variance at each level
  Q<- LKrig.precision( LKinfo)
  Q<- as.matrix( Q)
  look<- (X)%*% solve( Q) %*%t(X)
  marginal.var<- sum(unlist(LKinfo$alpha))
  test.for.zero( diag(look), rep( marginal.var,10),
                tag="normalization to unit variance at each level",
                tol=5e-7)
  look2<- LKrig.cov( x1, LKinfo= LKinfo, marginal=TRUE)
  test.for.zero( look2, rep(marginal.var,10) ,
                tag="normalization based on logic in LKrig.cov",
                tol=5e-7)
# check full covariance matrix
  look3<- LKrig.cov( x1, x1,LKinfo= LKinfo)
  test.for.zero( look3, look, tag="full covariance from matrix expressions")
# Now test w/o normalization
  LKinfo$normalize<- FALSE
  X<- LKrig.basis(x1, LKinfo)
  X<- as.matrix(X)
# check on normalization to unit variance at each level
  Q<- LKrig.precision( LKinfo)
  Q<- as.matrix( Q)
  look<- (X)%*% solve( Q) %*%t(X)
  look3<- LKrig.cov( x1, x1,LKinfo= LKinfo)
  test.for.zero( look3, look,
             tag="full covariance from matrix expressions w/o normalization")
#
# check of component covariance matrices
   LKinfo<- LKrigSetup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=5,
                        a.wght=c(5,6,7), alpha=c(4,2,1), NC.buffer=0)
  set.seed(123)
  x1<- cbind( runif( 10), runif(10))
  x2<- cbind(0,0)
  comp<- matrix( NA,10, 3)   
  for ( l in 1:3){
       LKinfo.temp<- LKrigSetup( cbind( c( -1,1), c( -1,1) ), nlevel=1, NC=LKinfo$latticeInfo$mx[l],
                        a.wght=c(5,6,7)[l], alpha=1.0, NC.buffer=0)
       comp[,l]<- LKrig.cov(x1,x2,LKinfo.temp )
  }
  look1<- comp%*%c(unlist( LKinfo$alpha))
  look3<- LKrig.cov( x1, x2, LKinfo= LKinfo)
  test.for.zero( look1, look3, tag="comp normalized cov and LKrig.cov")
#
# check construction with spatial a.wght
  cat("Now check spatial a.wght and alpha", fill=TRUE)

  LKinfo<- LKrigSetup( cbind( c( -1,1), c( -1,1) ), nlevel=1, NC=5,
                        a.wght=4,
                        alpha=1)
  a.wght<-  list( matrix(4 + (1:LKinfo$latticeInfo$m)*.1, LKinfo$latticeInfo$mx[1,2], LKinfo$latticeInfo$mx[1,2]))
  LKinfo<- LKrigSetup( cbind( c( -1,1), c( -1,1) ), nlevel=1, NC=5,
                        a.wght=a.wght,
                        alpha=1)

  look<- LKrig.precision( LKinfo=LKinfo, return.B=TRUE)
  look2<- as.matrix( look)
  test.for.zero( diag( look2), a.wght[[1]], tag="spatial a.wght 1 level")
# three levels
  LKinfo0 <- LKrigSetup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=4,
                        a.wght=c(5,5,5),
                        alpha=c(1,1,1), edge=FALSE)
  N<- LKinfo0$latticeInfo$mx[,1]*LKinfo0$latticeInfo$mx[,2]
  a.wght<-  list(
                 matrix(4 +  (1:N[1])*.1, LKinfo0$latticeInfo$mx[1,1], LKinfo0$latticeInfo$mx[1,2] ),
                 matrix(4 +  (1:N[2])*.1, LKinfo0$latticeInfo$mx[2,1], LKinfo0$latticeInfo$mx[2,2] ),
                 matrix(4 +  (1:N[3])*.1, LKinfo0$latticeInfo$mx[3,1], LKinfo0$latticeInfo$mx[3,2] )
                )
   LKinfo <- LKrigSetup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=4,
                        a.wght=a.wght,
                        alpha=c(1,1,1), edge=FALSE)
  look<- LKrig.precision( LKinfo=LKinfo, return.B=TRUE)
  look2<- as.matrix( look)
  test.for.zero( diag( look2), unlist(a.wght), tag="spatial a.wght 3 levels")
#
# STOPPPED HERE

# three levels nonzero buffer
  LKinfo0 <- LKrigSetup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=4,NC.buffer=3,
                        a.wght=c(5,5,5),
                        alpha=c(1,1,1), edge=FALSE)
  N<- LKinfo0$latticeInfo$mx[,1]*LKinfo0$latticeInfo$mx[,2]
  alpha<-  list(  (1:N[1])*.1, (1:N[2])*.1, (1:N[3])*.1)
  LKinfo <- LKrigSetup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=4, NC.buffer=3,
                        a.wght=5,
                        alpha=alpha, edge=FALSE)
  look<- LKrig.precision( LKinfo, return.B=TRUE)
  look2<- as.matrix( look)
 
  LKinfo2 <- LKrigSetup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=4, NC.buffer=3,
                        a.wght=5,
                        alpha=c(1,1,1), edge=FALSE)
  look3<- as.matrix(LKrig.precision( LKinfo2, return.B=TRUE))
  look3<-  diag( 1/sqrt(unlist(alpha)))%*%look3
  test.for.zero( look3, look2, tag=" 3 levels spatial alpha buffer=0")
  look4<-  as.matrix(LKrig.precision( LKinfo))
  test.for.zero( t(look3)%*%look3, look4, tag="3 levels spatial alpha Q buffer=3")


# test of buffer grid points.
 NC.buffer<- 7
  LKinfo0 <- LKrigSetup( cbind( c( -1,1), c( 0,3) ), nlevel=3, NC=12,
                        a.wght=c(5,6,7), alpha=c(4,2,1), NC.buffer=0)
  LKinfo7 <- LKrigSetup( cbind( c( -1,1), c( 0,3) ), nlevel=3, NC=12,
                        a.wght=c(5,6,7), alpha=c(4,2,1), NC.buffer=NC.buffer)
  check.sum <- 0
  for(k in 1:3){
    check.sum <- check.sum + length((LKinfo0$latticeInfo$grid[[k]])$x) +2*NC.buffer - length((LKinfo7$latticeInfo$grid[[k]])$x)
    check.sum <- check.sum + length((LKinfo0$latticeInfo$grid[[k]])$y) +2*NC.buffer - length((LKinfo7$latticeInfo$grid[[k]])$y)
    
  }
  test.for.zero( check.sum, 0, relative=FALSE, tag="adding  7 buffer points")
# getting the margins and spatial domain right
# this test works because x y aspects are both divided evenly by delta.
   
 check.sum<-0
 for(k in 1:3){
    check.sum <- check.sum + (LKinfo0$latticeInfo$grid[[k]])$x[1] - (LKinfo7$latticeInfo$grid[[k]])$x[NC.buffer +1] 
    check.sum <- check.sum +  (LKinfo0$latticeInfo$grid[[k]])$y[1] - (LKinfo7$latticeInfo$grid[[k]])$y[NC.buffer +1]
    m2<- LKinfo7$mx[k,1]
    n2<- LKinfo7$mx[k,2]
    check.sum <- check.sum + max((LKinfo0$latticeInfo$grid[[k]])$x) - (LKinfo7$latticeInfo$grid[[k]])$x[m2 - NC.buffer] 
    check.sum <- check.sum + max((LKinfo0$latticeInfo$grid[[k]])$y) - (LKinfo7$latticeInfo$grid[[k]])$y[n2 -NC.buffer] 
  }
 test.for.zero( check.sum, 0, relative=FALSE,tol=1e-8, tag="correct nesting of buffer=0  grid")
 cat("Done testing LKrig.precision",fill=TRUE)
  options( echo=FALSE)












