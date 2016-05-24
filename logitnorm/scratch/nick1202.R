# email of Jarle Tufo

momentsLogitnorm(5,0.1)
#        mean          var 
#1.950053e-11 1.943755e-11 

#mean          var 
#9.932743e-01 4.484069e-07		# with new default abs.tol=0 


momentsLogitnorm(5,0.1,subdivisions=10000, rel.tol=max(500*.Machine$double.eps, 0.5e-28), abs.tol=0)
momentsLogitnorm(5,0.1,subdivisions=10000, abs.tol=0)


library(lhs)
qfunct <- function(mean,sd) {    
	zz <- qnorm(randomLHS(100000,1),mean,sd)
	zz <- exp(zz)/(1+exp(zz))
	c(mean(zz),var(zz)) }

qfunct(5,0.1)
#[1] 9.933e-01 4.484e-07


tmp <- plogis(rnorm(1000000,5,0.1))
c( mean(tmp), var(tmp) )
# 0.993273		# confirmes the qfunct solution

mean(plogis(rnorm(1000000,5,0.1))) 

qfunct(4,1)
