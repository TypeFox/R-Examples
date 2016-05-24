# ex2.16.R
set.seed(123)
d <- expression(-4*tanh(2*x))
d.x <- expression(-(4 * (2/cosh(2 * x)^2)))
A <- function(x) -(0.5+6/4)*log(cosh(2*x))
X0 <- rt(1, df=4)/2
F <- function(x) log(x + sqrt(1+x^2))/2
Y0 <- F(X0) 
sde.sim(method="EA", delta=1/20, X0=Y0, N=500, drift=d, 
      drift.x=d.x, A=A, k1=-4,k2=8) -> Y
X <- sinh(Y)
ts(cbind(X,Y),start=0,delta=1/20) -> XY
plot(XY,main="Original scale X vs transformed Y")
