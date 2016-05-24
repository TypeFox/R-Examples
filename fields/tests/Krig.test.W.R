# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library(fields)
options( echo=FALSE)
test.for.zero.flag<- 1
#
#
#  test of off diagonal weight matrix for obs
#  Check against linear algebra
#
#cat("A very nasty case with off diagonal weights",fill=TRUE)

set.seed(123)
x<- matrix( runif( 30), 15,2)
y<- rnorm( 15)*.01 + x[,1]**2 +  x[,2]**2

#weights<- rep( 1, 15)

weights<- runif(15)*10


# WBW
# double check that just diagonals work. 

lambda.test<- .6
Krig( x,y,cov.function=Exp.cov,weights=weights)-> out
Krig( x,y,cov.function=Exp.cov,weights=weights, lambda=lambda.test)-> out2
Krig.coef( out, lambda=lambda.test)-> test

W<- diag( weights)
W2<- diag( sqrt(weights))


K<- Exp.cov(x,x) 
M<- (lambda.test*solve(W)  + K);T<- fields.mkpoly(x, 2)
temp.d<- c(solve( t(T) %*% solve( M)%*%T) %*% t(T)%*% solve( M) %*%y)
temp.c<- solve( M)%*% (y - T%*% temp.d)
#

# test for d coefficients
test.for.zero( test$d, out2$d, tag=" d coef diag W fixed lam")
test.for.zero( temp.d, out2$d, tag=" d coef diag W")
# test for c coefficents
test.for.zero( test$c, out2$c, tag="c coef diag W fixed lam" )
test.for.zero( temp.c, out2$c, tag="c coef  diag W " )



# the full monty

temp.wght<- function(x, alpha=.1){
  Exp.cov( x, theta=alpha) }

Krig( x,y,
     cov.function=Exp.cov,weights=weights, wght.function= temp.wght,
    )-> out.new

W2<-out.new$W2
W<- out.new$W



test.for.zero( c( W2%*%W2), c( W), tag=" sqrt of W")

Krig( x,y, cov.function=Exp.cov,weights=weights, W= out.new$W)-> temp

test.for.zero( predict(temp, y= y), predict(out.new,y=y), 
tag=" Test of passing W explicitly")



K<- Exp.cov(x,x); lambda.test<- .6; 
M<- (lambda.test*solve(W)  + K);T<- fields.mkpoly(x, 2)
temp.d<- c(solve( t(T) %*% solve( M)%*%T) %*% t(T)%*% solve( M) %*%y)
temp.c<- solve( M)%*% (y - T%*% temp.d)
# 
Krig.coef( out.new,lambda=lambda.test )-> out2

test.for.zero( temp.d, out2$d, tag=" d coef full W")
# test for c coefficents
test.for.zero( temp.c, out2$c, tag="c coef full W" )


####
### testing the GCV function 

lambda<- out.new$lambda

Krig.Amatrix( out.new, lambda=lambda)-> Alam

test.for.zero( Alam%*%y , predict(out.new), tag="A matrix")

N<- length( y)
test<- sum( diag( Alam))
# compare to
test2<- out.new$eff.df

test.for.zero( test,test2, tag=" check trace of A")

Krig.fgcv.one( lam=lambda, out.new)-> test
# compare to
test2<- (1/N)*sum(  
               (out.new$W2%*%(y - c(Alam%*% y) ))**2 
                               ) / (1- sum(diag( Alam))/N)**2

test.for.zero( test,test2,tol=.5e-7, tag="GCV one" )

cat( "all done  testing off diag W case", fill=TRUE)
options( echo=TRUE)
