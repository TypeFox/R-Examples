if(!require("GNE"))stop("this test requires package GNE.")

limphiFB <- function(x, ab) 
	( sqrt(x^2*ab[1]^2 + x^2*ab[2]^2) - x*(ab[1]+ab[2]) )/x
	
	
x <- seq(-.1, .1, .001)
n <- 10

plot(x, limphiFB(x, c(1,1)), type="l", ylim=c(-5, 2), 
	main="limits of phiFB(ta, tb)/t as t-> 0", xlab="t", ylab="phiFB(ta, tb)/t")

for(i in 1:n)
lines(x, limphiFB(x, c(1/i,1/(i-1))), col=grey( (i+5)/(n+15) ) )

for(i in 1:n)
lines(x, limphiFB(x, c(0,1/(i-1))), col=grey( (i+5)/(n+5) ) )
	
	
for(i in 1:n)
lines(x, limphiFB(x, c(1,1/(i-1))), col=grey( (i+5)/(n+5) ) )
	


somevals <- sapply(1:n, function(i) limphiFB(x, c(3,1/(i+1))))	
apply(somevals, 2, range, na.rm=TRUE)
	
	
	

limphigraFB	<- function(x, ab)
	sign(x)*ab[1]/sqrt(sum(ab^2)) - x*sum(ab)
	
		
x <- seq(.001, .1, .001)
n <- 10

plot(x, limphigraFB(x, c(1,1)), type="l", ylim=c(-1, 1), xlim=range(-x, x), 
	main="limits of GrAphiFB(ta, tb) as t-> 0", xlab="t", ylab="GrAphiFB(ta, tb)")
lines(-x, limphigraFB(-x, c(1,1)))

for(i in 1:n)
{
lines(x, limphigraFB(x, c(1/i,1/(i-1))), col=grey( (i+5)/(n+15) ) )
lines(-x, limphigraFB(-x, c(1/i,1/(i-1))), col=grey( (i+5)/(n+15) ) )
}


for(i in 1:n)
{
lines(x, limphigraFB(x, c(1,1/(i-1))), col=grey( (i+5)/(n+15) ) )
lines(-x, limphigraFB(-x, c(1,1/(i-1))), col=grey( (i+5)/(n+15) ) )
}
		
for(i in 1:n)
{
lines(x, limphigraFB(x, c(10,1/(i-1))), col=grey( (i+5)/(n+15) ) )
lines(-x, limphigraFB(-x, c(10,1/(i-1))), col=grey( (i+5)/(n+15) ) )
}