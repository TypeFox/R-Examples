# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library(fields)
#
#
#  test of fixed lambda case
#  Check against linear algebra
#

options( echo=FALSE)
test.for.zero.flag<-1

Krig( ChicagoO3$x, ChicagoO3$y, theta=50)-> fit

x<- ChicagoO3$x
K<- Exp.cov(x, x,theta=50)
T<- fields.mkpoly(x, 2)
W<- diag( 20)
 lambda<- fit$lambda
M<- (lambda* diag(20) + K) 
###########################
test.d<- c(solve( t(T) %*% solve( M)%*%T) %*% t(T)%*% solve( M) %*% fit$yM)
test.c<- solve( M)%*% ( fit$yM - T%*% test.d)

#compare to  fit$d
test.for.zero( test.d, fit$d, tag="Compare d coef" )
#compare to  fit$d
test.for.zero( test.c, fit$c,tag="Compare c coef" )

Krig( ChicagoO3$x, ChicagoO3$y, theta=50,lambda= fit$lambda)-> fit2
#compare to  fit$d
test.for.zero( test.d, fit2$d, tag="Compare d coef fixed lambda" )
#compare to  fit$d
test.for.zero( test.c, fit2$c,tag="Compare c coef fixed lambda" )

# test of Krig.coef

Krig.coef( fit)->test
test.for.zero( test.d, test$d, tag="d coef Krig.coef" )
test.for.zero( test.c, test$c, tag= "c coef Krig.coef" )

Krig.coef( fit2)->test
test.for.zero( test.d, test$d,tag="d coef Krig.coef fixed" )
test.for.zero( test.c, test$c, tag="c coef Krig.coef fixed" )
# checking A matrix in the case of noreps

set.seed( 222)
weights<-  10+ runif( length(ChicagoO3$y))
#weights<- rep( 1, 20)
test2<- Krig( ChicagoO3$x, ChicagoO3$y, theta=50, weights= weights)
Atest<- Krig.Amatrix( test2)
K<-Exp.cov(ChicagoO3$x, ChicagoO3$x,theta=50)
H<- matrix(0, 23,23)
H[(1:20)+3 , (1:20)+3]<- K
X<- cbind( fields.mkpoly( ChicagoO3$x, 2), K)
lambda<- test2$lambda
 Alam <-  X%*%solve(
                 t(X)%*%diag(weights)%*%X + lambda*H
                 )%*% t(X)%*%diag(weights) 
 test.for.zero( Alam, Atest, tag="Amatrix no reps", tol=5e-8)

# test for new y fixed case
set.seed( 123)
ynew<- rnorm( fit2$N)

test.d<- c(solve( t(T) %*% solve( M)%*%T) %*% t(T)%*% solve( M) %*% ynew)
test.c<- solve( M)%*% ( ynew - T%*% test.d)

Krig.coef( fit, y= ynew)->test
test.for.zero( test.d, test$d, tag= "d coef new y" )
test.for.zero( test.c, test$c, tag="c coef new y" )


Krig.coef( fit2, y= ynew)->test
test.for.zero( test.d, test$d, tag= "d coef new y fixed" )
test.for.zero( test.c, test$c, tag=" c coef new y fixed"  )

# test for multiple new y's
Krig.coef( fit2, y= cbind( ynew+ rnorm(fit2$N), ynew))->test2
test.for.zero( test.d, test2$d[,2], tag= "d coef several new y fixed" )
test.for.zero( test.c, test2$c[,2], tag=" c coef several new y fixed"  )


#cat("done with simple Krig data", fill=TRUE)


# These tests are about whether decompositions 
# handle just a fixed lambda or are more general 

# checking passing lambda or df to Krig

Tps( ChicagoO3$x, ChicagoO3$y,lambda=.001 )-> out
predict( out, lambda=.001)-> out2
test.for.zero( out2, predict( out), tag="Tps with fixed lam")

Tps( ChicagoO3$x, ChicagoO3$y, df=5)-> out
predict( out, df=5)-> out2
test.for.zero( out2, predict( out), tag="Tps with fixed df")

# same for Krig

Krig( ChicagoO3$x, ChicagoO3$y, theta=50,lambda=.5)-> out0
Krig( ChicagoO3$x, ChicagoO3$y, theta=50,lambda=.5,GCV=TRUE)-> out
test.for.zero( 
      predict(out0), predict( out), tag="Krig with fixed lam argument")

Krig( ChicagoO3$x, ChicagoO3$y, theta=50)-> out0
Krig( ChicagoO3$x, ChicagoO3$y, theta=50, df=6,GCV=TRUE)-> out
predict( out0, df=6)-> out2
test.for.zero( out2, predict( out), tag="Krig with fixed lam argument")


#cat("A very nasty case with knots and weights",fill=TRUE)

set.seed(123)
x<- matrix( runif( 30), 15,2)
y<- rnorm( 15)*.01 + x[,1]**2 +  x[,2]**2
knots<- x[1:5,]
weights<- runif(15)*10

# compare to 
Krig( x,y, knots=knots, cov.function=Exp.cov,weights=weights)-> out.new
Krig( x,y, knots=knots, cov.function=Exp.cov,weights=weights, 
          lambda=1)-> out.new2

# compute test using linear algebra

K<- Exp.cov( knots, knots)
H<- matrix(0, 8,8)
H[4:8, 4:8]<- K
X<- cbind( fields.mkpoly( x, 2), Exp.cov( x, knots))
lambda<-1


c(   solve(t(X)%*%(weights*X) + lambda*H)%*% t(X)%*% (weights*y) )-> temp
temp.c<- temp[4:8]
temp.d<- temp[1:3]


# test for d coefficients
test.for.zero( out.new2$d, temp.d, tag=" d coef")
# test for c coefficents
test.for.zero( out.new2$c, temp.c, tag="c coef" )


# compare to 
Krig.coef( out.new, lambda=1)->test
# and


# test for d coefficients
test.for.zero( temp.d, test$d, tag="d new y Krig.coef")
# test for c coefficents
test.for.zero( temp.c, test$c, tag="c new y Krig.coef" )


# and 
Krig.coef( out.new2, lambda=1)-> test

# test for d coefficients
test.for.zero( temp.d, test$d, tag= "d fixed case")
# test for c coefficents 
test.for.zero( temp.c, test$c, tag=" c fixed case" )



#cat( "done with knots and weights case", fill=TRUE)

#
#  test with new y
#  

lam.test <- 1.0

ynew<- 1:15

Krig( x,y, knots=knots, cov.function=Exp.cov,weights=weights)-> out.new
Krig( x,y, knots=knots, cov.function=Exp.cov,weights=weights, 
                 lambda=lam.test)-> out.new2
### compare to 
##Krig( x,ynew, knots=knots, cov.function=Exp.cov,weights=weights)-> out.new
##Krig( x,ynew, knots=knots, cov.function=Exp.cov,weights=weights, 
##                 lambda=lam.test)-> out.new2

c(   solve(t(X)%*%(weights*X) + lam.test*H)%*% t(X)%*% (weights*ynew) )-> temp
temp.d<- temp[1:3]
temp.c<- temp[4:8]

#compare 
Krig.coef( out.new,lambda=lam.test,y=ynew)-> test

# test for d coefficients
test.for.zero( temp.d, test$d, tag=" d new y")
# test for c coefficents 
test.for.zero( temp.c, test$c,tag= "c new y" )


Krig.coef( out.new2,y=ynew)-> test

# test for d coefficients
test.for.zero( temp.d, test$d, tag= "d new y fixed")
# test for c coefficents 
test.for.zero( temp.c, test$c, tag= "c new y fixed" )



#cat( "done with new y case for nasty data ", fill=TRUE)


#
#cat("test with reps" , fill=TRUE)
#

set.seed(133)
x<- matrix( runif( 30), 15,2)*2
x<- rbind( x,x, x[3:7,])
y<- rnorm( nrow( x))*.05 + + x[,1]**2 +  x[,2]**2
# perturb so that this example does not generate (harmless) warnings in gcv search
y[20] <- y[20] + 1
weights<- runif( nrow( x))*10 
knots<- x[1:10,]

Krig( x,y, knots=knots,  weights=weights, cov.function=Exp.cov)-> out.new



lambda<- 1.0
NP<- out.new$np
NK <- nrow( knots)
K<- Exp.cov( knots, knots)
H<- matrix(0, NP,NP)
H[(1:NK)+3 , (1:NK)+3]<- K
X<- cbind( fields.mkpoly( x, 2), Exp.cov( x, knots))

# compare to 
test<- c(   solve(t(X)%*%diag(weights)%*%X + lambda*H)%*% 
t(X)%*%diag(weights)%*% y )

test[1:3]-> temp.d
test[(1:NK)+3]-> temp.c

Krig( x,y, knots=knots,  weights=weights,lambda=lambda,
 cov.function=Exp.cov)-> out.new

# test for d coefficients
test.for.zero( temp.d, out.new$d, tag=" d reps")
# test for c coefficents 
test.for.zero( temp.c, out.new$c, tag="c reps" )


Krig( x,y, knots=knots,  weights=weights, cov.function=Exp.cov)-> out.new

#compare to
test<-  sum(weights*
     (y-X%*%solve(t(X)%*%diag(weights)%*%X) %*% t(X)%*%diag(weights)%*% y)**2
    )

test.for.zero(out.new$pure.ss, test, tag=" pure sums of squares")



#cat("done with reps case", fill=TRUE)

##################################
#cat( "test  A matrix",fill=TRUE)
##################################
 
set.seed(133)
x<- matrix( runif( 30), 15,2)*2  
x<- rbind( x,x, x[3:7,])
y<- rnorm( nrow( x))*.05 + + x[,1]**2 +  x[,2]**2
# perturb so that this example does not generate (harmless) warnings in gcv search
y[20] <- y[20] + 1
weights<- runif( nrow( x))*10
knots<- x[1:10,]

Krig( x,y, knots=knots,  weights=weights, cov.function=Exp.cov)-> out.new

NP<- out.new$np
NK <- nrow( knots)
K<- Exp.cov( knots, knots)
H<- matrix(0, NP,NP)
H[(1:NK)+3 , (1:NK)+3]<- K
X<- cbind( fields.mkpoly( x, 2), Exp.cov( x, knots))



lambda<- out.new$lambda
 Alam= X%*%solve(t(X)%*%diag(weights)%*%X + lambda*H)%*% t(X)%*%diag(weights)
 
test<- c(Alam%*% y)
# compare to
test2<-predict( out.new)

test.for.zero( test,test2, tag="Amatrix prediction")

#
test<- sum( diag( Alam))
test2<- out.new$eff.df
     
test.for.zero( test,test2)

Krig.Amatrix( out.new, lambda=lambda)-> Atest
test.for.zero( sum( diag(Atest)),test2, tag=" trace from A matrix")

test.for.zero( Atest%*%out.new$yM, predict(out.new))

yjunk<- rnorm( 35)
yMtemp<- Krig.ynew(out.new, yjunk)$yM
test.for.zero( Atest%*%yMtemp, predict(out.new, y=yjunk),
tag="A matrix predict with new y")

test.for.zero( Atest%*%yMtemp, predict(out.new, yM= yMtemp), 
tag="A matrix predict compared to collapsed yM")


test.pure.ss<-  sum(weights*
     (y-X%*%solve(t(X)%*%diag(weights)%*%X) %*% t(X)%*%diag(weights)%*% y)**2
    ) 


test.for.zero( out.new$pure.ss, test.pure.ss,tag="pure sums of squares")

#cat("done with A matrix case", fill=TRUE)
#
# check of GCV etc. 

lambda<- out.new$lambda
 Alam= X%*%solve(t(X)%*%diag(weights)%*%X + lambda*H)%*% t(X)%*%diag(weights)

test<- c(Alam%*% y)
# compare to 
test2<-predict( out.new)

#test.for.zero( test,test2, tag="double check A matrix predict")


N<- length( y)
test<- sum( diag( Alam))
# compare to 
test2<- out.new$eff.df

test.for.zero( test,test2, tag=" check trace")

