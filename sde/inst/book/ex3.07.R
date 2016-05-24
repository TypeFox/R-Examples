# require(sde)

# ex3.07.R
d1 <- function(t,x,theta) theta[1]*(theta[2]-x)
s1 <- function(t,x,theta) theta[3]*sqrt(x)
from <- 0.08
x <- seq(0,0.2, length=100)
sle10 <- NULL
sle2 <- NULL
sle5 <- NULL
true <- NULL
set.seed(123)
for(to in x){
  sle2 <- c(sle2, dcSim(from, to, 0.5, d1, s1, theta=c(2,0.02,0.15), 
                  M=50000,N=2))
  sle5 <- c(sle5, dcSim(from, to, 0.5, d1, s1, theta=c(2,0.02,0.15), 
                  M=50000,N=5))
  sle10 <- c(sle10, dcSim(from, to, 0.5, d1, s1, theta=c(2,0.02,0.15), 
                  M=50000,N=10))
  true <- c(true, dcCIR(to, 0.5, from, c(2*0.02,2,0.15)))
}
par(mar=c(5,5,1,1))
plot(x, true, type="l", ylab="conditional density")
lines(x, sle2, lty=4)
lines(x, sle5, lty=2)
lines(x, sle10, lty=3)
legend(0.15,20, legend=c("exact","N=2", "N=5", "N=10"), lty=c(1,2,4,3))
