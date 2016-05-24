# empirical and theoretical normal cumulative plots

cdf.plot <- function(x){
 # number of observations
 nx <- length(x)
 # data sorted	
 vx <- sort(x)
 # ranks as fractions (quantiles)
 Fx <- seq(nx)/nx
 # empirical cumulative plot	
 plot(vx, Fx, xlab="x", ylab="F(x)")
 # mean and stdev 	
 mx <- mean(vx); sdx <- sd(vx)
 # do sequence with increased res 	
 vnx <- seq(vx[1],vx[nx],length=100)
 # theoretical cdf values 
 Fxn <- pnorm(vnx,mx,sdx)
 # theoretical cumulative plot CDF 	
 lines(vnx,Fxn) 		
 ret.val <- list(Fx.Emp=Fx,Fx.Theo=Fxn)
}

# Exploratory Data Analysis (EDA)

eda6 <- function(x,label="x"){
 # arguments x = array, label= string to label axis 
 # divide graphics window in 6 panels
 panel6(size=7)
 # index plot
 plot(x, ylab=label, main="Index plot",cex.main=0.7) 
 # box-whisker plot	
 boxplot(x,ylab=label); title("Boxplot",,cex.main=0.7)	
 # histogram
 hist(x, main="Histogram",xlab=label,cex.main=0.7)	
 #density approx to histogram
 plot(density(x), main="Density approximation",xlab=label,cex.main=0.7)
 # quantile-quantile plot require package car
 qqPlot(x,col.lines=1,lwd=1,grid=F)
 title("QQ plot",cex.main=0.7)
 # empirical cumulative plot
 # standardize observations before plot
 z <- (x-mean(x))/sqrt(var(x))
 plot(ecdf(z), xlab=paste("Standardized ", label,sep=""),
               main= "ECDF vs Std. normal",cex.main=0.7)
 # generate standard normal and plot 
 Z <- seq(-4,+4,0.1); lines(Z, pnorm(Z,0,1))
}

# panel functions

panel2 <- function (size,int="r"){
mat <- matrix(1:2,2,1,byrow=T)
layout(mat, widths=rep(size,2), heights=rep(size/2,2), TRUE)
par(mar=c(4,4,1,.5),xaxs=int, yaxs=int)
}
panel4 <- function (size,int="r"){
mat <- matrix(1:4,2,2,byrow=T)
layout(mat, widths=rep(size/2,2), heights=rep(size/2,2), TRUE)
par(mar=c(4,4,1,.5),xaxs=int, yaxs=int)
}
panel6 <- function (size,int="r"){
mat <- matrix(1:6,3,2,byrow=T)
layout(mat, widths=rep(size/2,2), heights=rep(size/3,3), TRUE)
par(mar=c(4,4,1,.5),xaxs=int, yaxs=int)
}

panel3 <- function (size,int="r"){
mat <- matrix(1:3,3,1,byrow=T)
layout(mat, widths=rep(size,3), heights=rep(size/3,3), TRUE)
par(mar=c(4,4,1,.5),xaxs=int, yaxs=int)
}

 



