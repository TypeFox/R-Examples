library( LatticeKrig)
#
options(echo=FALSE)
test.for.zero.flag<- 1
#
# Near interpolation of a 3d function 
set.seed( 123)
N<- 3e4
x<-  matrix( runif(3* N,-1,1), ncol=3, nrow=N)
y<-   10*exp( -rdist( x, rbind( c(.5,.5,.6) ) )/.5)
glist<- list( x1=seq( -1,1,,30), x2=seq( -1,1,,30), x3= 0)
xgrid<- make.surface.grid( glist)
yTrue<-   10*exp( -rdist( xgrid, rbind( c(.5,.5,.6) ) )/.5)

LKinfo<- LKrigSetup( x,  nlevel=1,  a.wght= 6.2, NC=8, NC.buffer=2,
                    LKGeometry="LKBox", normalize=TRUE, mean.neighbor=200,
                    choleskyMemory=list(nnzR= 2e6))

out1<- LatticeKrig( x,y, LKinfo=LKinfo)
yTest<- predict( out1, xgrid)
# accuracy within a few percent relative error.
test.for.zero( mean( abs(yTest- yTrue)/yTrue), 0, relative=FALSE, tol=3e-2 )

# Larger model

LKinfo<- LKrigSetup( x,  nlevel=3,  a.wght= 8, alpha=c( 1,.5, .2),
                    NC=3, NC.buffer=1,
                    LKGeometry="LKBox", normalize=TRUE, mean.neighbor=200,
                    choleskyMemory=list(nnzR= 2e6))

# test of finding nearest lattice neighbors
remove( test.for.zero.flag)
   cat(" Exahautive test of 3d box lattice", fill=TRUE)
for(  level in 1:LKinfo$nlevel){
m1<- LKinfo$latticeInfo$mLevel[level]
mx1<- LKinfo$latticeInfo$mx[level,]
indexgrid<- as.matrix( expand.grid( list(1:mx1[1], 1:mx1[2], 1:mx1[3] )))
look<- LKrigSAR( LKinfo, Level=level)
BigD<- rdist( indexgrid, indexgrid)
allNodes<- 1:m1
for(  k in 1:m1){
#	cat( k, " ")
	i1<- look$ind[ look$ind[,1]==k,2]
	i2<- allNodes[ BigD[k,]<=1.0]
#	print( length( i1) - length(i2))
	test.for.zero( sort( i1), sort( i2)) 	
 }
    cat(" ", fill=TRUE)
    cat( " done  with level ", level, fill=TRUE)
}



glist<- list( x1=seq( -1,1,,30), x2=seq( -1,1,,30), x3= 0)
xgrid<- make.surface.grid( glist)

look3<- LKrig.cov( xgrid, LKinfo=LKinfo, marginal=TRUE)
varTest<- sum( unlist(LKinfo$alpha))
test.for.zero( look3, varTest, tag="checking computation of marginal variance LKBox") 
 






