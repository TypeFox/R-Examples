source("../scripts/RK4.r")

LV <- function(t = NULL, x, a, b, g, e, K=Inf)
  c( a*x[1]*(1 - x[1]/K) - b*x[1]*x[2], g*b*x[1]*x[2] - e*x[2] )

h <- 0.005
x <- RK4(LV, 0, c(300,10), h, 500/h, a=0.05, b=0.0002, g=0.8, e=0.03, K=500)

#postscript(file="odes_LotkaVolterra.ps")
plot(x[,1], x[,2], type="l", xlab="prey", ylab="pred", main="Lotka-Volterra eqns")
#dev.off()
