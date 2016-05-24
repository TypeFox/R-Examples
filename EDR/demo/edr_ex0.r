require(EDR)
example1<-function(d,n,sigma){
f1 <- function(x) x^2
x <- matrix(2*rbeta(n*d,1,1)-1,n,d)
R1 <- matrix(0,d,d)
R1[1,] <- c(1,2,rep(0,d-2))/sqrt(5)
fx <- f1(x%*%t(R1)[,1])
y <- rnorm(n,fx,sigma)
list(y=y,x=x,fx=fx,R1=R1)
}
set.seed(1)
#
#  Estimate the effective dimension reduction space
#
d <- readline("Press 'Enter' for 25-dimensional example, otherwise provide the dimension:")

if(is.na(as.integer(d))) d <- 25 else d <- as.integer(d)
if(d<3) {
print("d should be larger than 2, d=25 is used \n")
d <- 25
}
if(d>50) { 
print("d  larger than 50 takes to much time for a demo, d=50 is used \n")
d <- 50
}
n <- 8*d
sigma <- .3
set.seed(1)
z<-example1(d,n,sigma)
zg<-list(x=z$x,y=z$y,fhat=z$fx,bhat=t(z$R1))
cat("Show projection into the effective dimension reduction space (EDR)\n")
plot.edr(zg,title="Projection into the true EDR",sm=FALSE)
readline("Press 'Enter return for random projections:")
oldpar<-par(mfrow=c(3,4),mar=c(1,1,3,.5),mgp=c(2,1,0))
for(i in 1:12){
zg$bhat <- svd(matrix(runif(d*d),d,d))$u
plot.edr(zg,title="Random projection",sm=FALSE)
}
readline("Press 'Enter return for estimation of the EDR:")
zedr <- edr(z$x,z$y,m=2,graph=TRUE,show=1)
readline("Press 'Enter return to summarize results:")

sedr <- summary(zedr,m=1,R=z$R1)
readline("Press 'Enter' to plot results:")
par(mfrow=c(1,2),mar=c(1,1,3,.5),mgp=c(2,1,0))
#zedr$fhat <- z$fx
plot(zedr,m=1,title="Projection into the estimated EDR")
zg$bhat<-t(z$R1)
zg$fhat <- z$fx
plot.edr(zg,title="Projection into the true EDR",sm=FALSE)
par(oldpar)
rm(z,zedr,sedr,d,n,sigma,example1)

