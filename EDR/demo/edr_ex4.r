require(EDR)
#
#        Example 2    (2D, d=10, n=100)
#
example2<-function(d,n,sigma){
f1b <- function(x) x[,1]*sin(sqrt(5)*x[,2])+x[,2]*sin(sqrt(5)*x[,1])
set.seed(1)
x <- matrix(2*rbeta(n*d,1,1)-1,n,d)
R1 <- matrix(0,d,d)
R1[1,] <- c(1,2,rep(0,d-2))/sqrt(5)
R1[2,] <- c(-2,1,2,rep(0,d-3))/3
fx <- f1b(x%*%t(R1)[,1:2])
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
if(d>50) { 
print("d  larger than 50 takes to much time for a demo, d=50 is used \n")
d <- 50
}
n <- 10*d
sigma <- 0.3
z<-example2(d,n,sigma)
cat("Run with graphical illustration (graph=TRUE)\n")
zedr <- edrcv(z$x,z$y,m=3,graph=TRUE,show=2,cvsize=min(10,d),hsm=seq(.1,.5,.01))
readline("Press 'Enter' to summarize results:")

sedr <- summary(zedr,m=2,R=z$R1)
readline("'Enter' to plot results:")

plot(zedr,m=2)
rm(z,zedr,sedr,d,n,sigma,example2)

