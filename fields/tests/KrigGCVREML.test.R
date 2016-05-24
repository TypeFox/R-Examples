library(fields)
#
#
#

options( echo=FALSE)
test.for.zero.flag<-1

############ various tests of GCV and REML
set.seed(133)
x0<- matrix( runif( 10*2), 10,2)*2  
x<- rbind( x0,x0, x0[3:7,])
y<- rnorm( nrow( x))*.05 + + x[,1]**2 +  x[,2]**2
weights<- 8 + runif( nrow( x))

# x0 are the unique values.


out.new<- Krig( x,y,   weights=weights, cov.function=Exp.cov)
n<- length(y)
n0<- nrow( x0)
NK <- nrow( x0) 
NP<- NK + 3
K<- Exp.cov( x0, x0)
H<- matrix(0, NP,NP)
H[(1:NK)+3 , (1:NK)+3]<- K
X<- cbind( fields.mkpoly( x, 2), Exp.cov( x, x0) )
X0<- cbind( fields.mkpoly( x0, 2), Exp.cov( x0, x0) )
Alam <-  X%*%solve(
                  t(X)%*%diag(weights)%*%X + out.new$lambda*H
                  )%*% t(X)%*%diag(weights) 
# predict sanity check using replicates 
set.seed( 123)
ynew<- rnorm(n)
test.for.zero( Alam%*%ynew, predict( out.new, y=ynew), tag=" predict sanity check",tol=3e-8) 

# predict using unique obs
ynew<- rnorm(nrow(x0))
Alam0<- X0%*%solve(
                  t(X0)%*%diag(out.new$weightsM)%*%X0 + out.new$lambda*H
                  )%*% t(X0)%*%diag(out.new$weightsM) 
 
# Alam0 is the A matrix                  
test.for.zero( Alam0%*%ynew, predict( out.new,x=x0, yM=ynew), tag="predict using direct linear algebra" )

#
test<- Krig.fgcv( lam=out.new$lambda, out.new)
y0<- out.new$yM
n0<- length(y0)
# compare to 
#test2<- (1/n0)*sum(  (y0 - c(Alam0%*% y0))**2 *out.new$weightsM) / (1- sum(diag( Alam0))/n0)**2
NUM<- mean(  (y0 - c(Alam0%*% y0))**2 *out.new$weightsM)  + out.new$pure.ss/( n -n0 )
DEN<- (1- sum(diag( Alam0))/n0)
test2<- NUM/ DEN^2
test.for.zero( test,test2, tag="GCV" )

test<- Krig.fgcv.one( lam=out.new$lambda, out.new)
N<- length(y) 
test2<- (1/N)*sum(  (y - c(Alam%*% y))**2 *weights) / 
                             (1- sum(diag( Alam))/N)**2                            
test.for.zero( test,test2, tag="GCV one" )

test<- Krig.fgcv.model( lam=out.new$lambda, out.new)  
y0<- out.new$yM
n0<- length(y0)
# compare to 
test2<- (1/n0)*sum(  (y0 - c(Alam0%*% y0))**2 *out.new$weightsM) / (1- sum(diag( Alam0))/n0)**2 + out.new$shat.pure.error**2
test.for.zero( test,test2,tag="GCV model")




####### tests with higher level gcv.Krig

data( ozone2)
x<- ozone2$lon.lat
y<- ozone2$y[16,]
Tps( x,y)-> out
gcv.Krig( out, tol=1e-10)-> out2

test.for.zero(out$lambda.est[1,-6], 
       out2$lambda.est[1,-6],tol=5e-4, tag="Tps/gcv for ozone2")

# try with "new" data (linear transform should give identical 
# results for GCV eff df

gcv.Krig( out, y=(11*out$y + 5), tol=1e-10 )-> out3

test.for.zero(out2$lambda.est[1,2], 
       out3$lambda.est[1,2],tol=1e-6, tag="Tps/gcv for ozone2 new data")

#cat("done with GCV case", fill=TRUE)



cat("done with GCV and REML tests", fill=TRUE)
options( echo=TRUE)