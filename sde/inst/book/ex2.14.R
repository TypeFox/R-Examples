# ex2.14.R
bY <- expression( (23-11*x^2+1.5*x^4-(x^6)/(2^4))/(2*x) ) 
bX <- expression( (6-11*x+6*x^2-x^3) ) 
sX <- expression( sqrt(x) )
 
set.seed(123)
X <- sde.sim(drift=bX, sigma=sX)
plot(X)
set.seed(123)
Y <- sde.sim(drift=bY, X0 = 2, method="shoji") 
plot((Y/2)^2)
