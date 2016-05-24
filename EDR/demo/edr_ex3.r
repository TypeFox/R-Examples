require(EDR)
#
#        Example 3  (2D, d=10, n=400)
#
example3<-function(d,n,sigma){
f1c <- function(x) sqrt(x[,1]^2+x[,2]^2)*sin(4*atan(x[,1]/(x[,2]+1e-5*sign(x[,2]))))
set.seed(1)
x <- matrix(2*rbeta(n*d,1,1)-1,n,d)
R1 <- matrix(0,d,d)
R1[1,] <- c(1,2,rep(0,d-2))/sqrt(5)
R1[2,] <- c(-2,1,2,rep(0,d-3))/3
fx <- f1c(x%*%t(R1)[,1:2])
y <- rnorm(n,fx,sigma)
list(y=y,x=x,fx=fx,R1=R1)
}
#
#  Estimate the effective dimension reduction space
#
d <- readline("Press 'Enter' for 10-dimensional example, otherwise provide the dimension:")

if(is.na(as.integer(d))) d <- 10 else d <- as.integer(d)
if(d<3) {
print("d should be larger than 2, d=10 is used \n")
d <- 10
}
if(d>25) { 
print("d  larger than 25 takes to much time for a demo, d=25 is used \n")
d <- 25
}
n <- 40*d
sigma <- 0.2
z<-example3(d,n,sigma)
cat("Run with additional diagnostics (graph=TRUE,trace=TRUE)\n")
zedr <- edr(z$x,z$y,m=3,graph=TRUE,show=2,trace=TRUE)
readline("Press 'Enter'  to summarize results:")

sedr <- summary(zedr,m=2,R=z$R1)
readline("Press 'Enter'  to plot results:")

plot(zedr,m=2)
readline("Press 'Enter'  to print results:")

print(zedr,m=2)
rm(z,zedr,sedr,d,n,sigma,example3)

