# LatticeKrig
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html


# test of radial basis function based on Wendland
# and using sparse formats
# Important check is of the FORTRAN function dfind2d
# that does pairwise distances among points within a specified range.

  library(LatticeKrig)
  options( echo=FALSE)
  test.for.zero.flag<-1

 
# check basic formula
 
  set.seed( 333)
  xLocation<- cbind( runif( 10, 3,5), runif( 10,3,5))
  xNew<- cbind( runif( 4, 3,5), runif( 4,3,5))
  LKinfo<- LKrigSetup( cbind( c(3,5), c(3,5)), NC=3, NC.buffer=2,
                     a.wght=5,  alpha=c(1), nlevel=1, normalize=TRUE)
  PHI<- LKrig.basis(xNew, LKinfo)
  Q<- LKrig.precision(LKinfo)
  covMatrix<- (PHI)%*% solve( Q)%*%t(PHI)
test.for.zero( diag(covMatrix), rep( 1, nrow( xNew)),
 tag="check using Qinverse formula", tol=1e-7)

covMatrix2<- LKrig.cov( xNew, LKinfo=LKinfo)                                       
test.for.zero( covMatrix,covMatrix2,
 tag="check using Qinverse formula full matrix",, tol=1e-7)

# multiple levels
 alpha<- c( 1, .8,.2)
 LKinfo<- LKrigSetup( cbind( c(3,5), c(3,5)), NC=2, NC.buffer=2,
                     a.wght=5,  alpha=alpha, nlevel=3, normalize=TRUE)
 varTest<- sum( alpha)                    
  PHI<- LKrig.basis(xNew, LKinfo)
  Q<- LKrig.precision(LKinfo)
  covMatrix<- (PHI)%*% solve( Q)%*%t(PHI)
test.for.zero( diag(covMatrix), rep( varTest, nrow( xNew)),
 tag=" Qinverse formula norm", tol=1e-7)

# multiple levels no normalization
 alpha<- c( 1, .8,.2)
 LKinfo<- LKrigSetup( cbind( c(3,5), c(3,5)), NC=2, NC.buffer=2,
                     a.wght=5,  alpha=alpha, nlevel=3, normalize=FALSE)                   
  PHI<- LKrig.basis(xNew, LKinfo)
  Q<- LKrig.precision(LKinfo)
  covMatrix<- (PHI)%*% solve( Q)%*%t(PHI)
  covMatrix2<- LKrig.cov( xNew, LKinfo=LKinfo)
  test.for.zero( covMatrix, covMatrix2,
               tag=" Qinverse formula nonorm")

# multiple levels
 alpha<- c( 1, .8,.2)
 LKinfo<- LKrigSetup( cbind( c(3,5), c(3,5)), NC=2, NC.buffer=2,
                     a.wght=8,  alpha=alpha, nlevel=3, normalize=TRUE,
                     BasisType="Tensor")
  varTest<- sum( alpha)                    
  PHI<- LKrig.basis(xNew, LKinfo)
  Q<- LKrig.precision(LKinfo)
  covMatrix<- (PHI)%*% solve( Q)%*%t(PHI)
  varTest<-LKrig.cov( xNew, marginal=TRUE, LKinfo=LKinfo)
  test.for.zero( diag(covMatrix), varTest,
                tag=" Qinverse formula Tensor no norm")
  test.for.zero( diag(covMatrix), sum(alpha),
                tag=" Qinverse formula Tensor no norm")


############################ end Q inverse formula checks

  set.seed( 333)
  xLocation<- cbind( runif( 10, 3,5), runif( 10,3,5))
  LKinfo<- LKrigSetup( cbind( c(3,5), c(3,5)), NC=3, NC.buffer=2,
                                        a.wght=5,  alpha=1, nlevel=1, normalize=TRUE)
                                                                                                                                                          
  LKinfo0<- LKrigSetup( cbind( c(3,5), c(3,5)), NC=3, NC.buffer=2,
                                        a.wght=5,  alpha=1, nlevel=1, normalize=FALSE)
  wght1<- LKrigNormalizeBasisFast.LKRectangle(LKinfo, Level=1, xLocation)
  wght0<-LKrig.cov( xLocation,   LKinfo= LKinfo0, marginal=TRUE)

  test.for.zero( wght0, wght1, tag=" 1 level Marginal variance compared to  fast normalize", tol=2e-7)
###  multiple levels

  alphaVec<- c( 1,.5,.2)
  LKinfo0<- LKrigSetup( cbind( c(3,5), c(3,5)), NC=3, NC.buffer=2,
                         a.wght=5, alpha=alphaVec, nlevel=3, normalize=FALSE )                         
  LKinfo<- LKrigSetup( cbind( c(3,5), c(3,5)), NC=3, NC.buffer=2,
                         a.wght=5, alpha=alphaVec,nlevel=3,normalize=TRUE)
                                        
  test1<-LKrigNormalizeBasisFast.LKRectangle( LKinfo, Level=1, xLocation )
  test2<-LKrigNormalizeBasisFast.LKRectangle( LKinfo, Level=2, xLocation )
  test3<-LKrigNormalizeBasisFast.LKRectangle( LKinfo, Level=3, xLocation )
 
  testVar1<- cbind( test1, test2, test3) %*% alphaVec
  testVar0<-LKrig.cov( xLocation,   LKinfo= LKinfo0, marginal=TRUE)
  test.for.zero( testVar0, testVar1, tag="Marginal variance and fast normalize", tol=1e-7)

cat( "Done with testing fast normalize algorithm", fill=TRUE)
options( echo=TRUE)

