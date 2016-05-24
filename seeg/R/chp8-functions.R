scan.geoeas.ppp <- function(filename){
# functions to scan geoeas files in point pattern
# miguel acevedo april 2004 
 n.coord.var <- scan(filename, skip=1, n=1)
 names.coord.var <- scan(filename, skip=2, what=(""), n=n.coord.var)
 coord.var <- read.table(filename, header=F, skip=n.coord.var+2, col.names=names.coord.var)
 return(coord.var)
}

# -------------------------
quad.chisq.ppp <- function(dataset,target.intensity){

#determine range
xr <- c( floor(min(dataset[,1])), ceiling(max(dataset[,1])) )
yr <- c( floor(min(dataset[,2])), ceiling(max(dataset[,2])) )

# make ppp
pppset <- ppp(dataset[,1],dataset[,2],xrange=xr,yrange=yr)

# describe it
summary(pppset)

# find n points
npoints <- length(dataset[,1])

# find n of quadrats 
n <- floor(sqrt(npoints/target.intensity))

# finding breakpoints
dxr <- (xr[2]-xr[1])/n
xbrks <- seq(xr[1],xr[2],dxr)
dyr <- (yr[2]-yr[1])/n
ybrks <- seq(yr[1],yr[2],dyr)
#count points in each cell with functions cut and table
 x <- cut(pppset$x, breaks=xbrks)
 y <- cut(pppset$y, breaks=ybrks)
 X <- table(y,x)

# Rows correspond to the y coordinate and the columns to the x coordinate
# because of the order we followed in the table command.

# Table X is converted to a regular matrix 
 Xm <- matrix(X,ncol=n,byrow=F)

# Note that larger y values will occur for the first rows.
# However, x and y have their origin in the SW corner,
# the y coordinate increases with decreasing row position and
# the x coordinate increases with column position.
# Thus, the sequence of the rows is such that the first row corresponds
# to the southernmost row and the last row corresponds to the northernmost row.
# The following inverts the order of the rows, so that the matrix overlays the map  
 Xint <- Xm; for(i in 1:n) Xint[n-(i-1),] <- Xm[i,]

 # estimate the min and max of the uniform; seq and expected prob
 # d <- sd(c(X))*sqrt(12)/2; a <- ceiling(mean(X)-d); b <- floor(mean(X)+d)
 # quantiles, prob and cumulative 
 #q <- seq(a, b, 1); dexp <- 1/(b-a); pexp <- seq(0,1,dexp)
 # odserved probs
 q <- seq(1,length(tabulate(X)),1)
 pobs <- tabulate(X)[q]/sum(tabulate(X))
 dexp <- 1/mean(X); pexp <- seq(0,1,1/(length(q)-1))


 # chisquare statistics
 X2 <- sum(((X-mean(X))^2)/ mean(X)); X2
 df <- n*n - 2; df
 p.value <- 1 - pchisq(X2, df)

  panel4(size=7);par(pty="s")
  plot(pppset$x,pppset$y, xlab="x",ylab="y")
  #draw grid lines 
  for(i in 1:(n+1))  abline(h=ybrks[i],lty=2)
  for(i in 1:(n+1))  abline(v=xbrks[i],lty=2)
  title("Point Pattern",cex.main=0.7)

 # use transpose of Xm. We could have applied the im(.) command.  
 image(xbrks, ybrks, t(Xm), col=gray((10:0)/10),xlab="x",ylab="y")
  for(i in 1:(n+1))  abline(h=ybrks[i])
  for(i in 1:(n+1))  abline(v=xbrks[i])
  title("Intensity",cex.main=0.7)
 # plot observed  
  plot(q,pobs, type="h",ylab="Proportion",xlab="Count",xlim=c(1,max(X)))
  title("Observed",cex.main=0.7)
  abline(h=dexp,lty=2,col=2)
  plot(ecdf(X), main="",ylab="F(Count)",xlab="Count",xlim=c(1,max(X))); lines(q,pexp,col=1)
  title("ECDF",cex.main=0.7)

return(list(pppset=pppset, Xint=Xint, intensity=mean(X), chisq=X2, p.value=p.value))
}

# -------------------------
quad.poisson.ppp <- function(dataset,target.intensity){

#determine range
xr <- c( floor(min(dataset[,1])), ceiling(max(dataset[,1])) )
yr <- c( floor(min(dataset[,2])), ceiling(max(dataset[,2])) )

# make ppp
pppset <- ppp(dataset[,1],dataset[,2],xrange=xr,yrange=yr)

# describe it
summary(pppset)

# find n points
npoints <- length(dataset[,1])

# find n of quadrats 
n <- floor(sqrt(npoints/target.intensity))

# finding breakpoints
dxr <- (xr[2]-xr[1])/n
xbrks <- seq(xr[1],xr[2],dxr)
dyr <- (yr[2]-yr[1])/n
ybrks <- seq(yr[1],yr[2],dyr)
#count points in each cell with functions cut and table
 x <- cut(pppset$x, breaks=xbrks)
 y <- cut(pppset$y, breaks=ybrks)
 X <- table(y,x)
# Rows correspond to the y coordinate and the columns to the x coordinate
# because of the order we followed in the table command.

# Table X is converted to a regular matrix 
 Xm <- matrix(X,ncol=n,byrow=F)

# Note that larger y values will occur for the first rows.
# However, x and y have their origin in the SW corner,
# the y coordinate increases with decreasing row position and
# the x coordinate increases with column position.
# Thus, the sequence of the rows is such that the first row corresponds
# to the southernmost row and the last row corresponds to the northernmost row.
# The following inverts the order of the rows, so that the matrix overlays the map  
 Xint <- Xm; for(i in 1:n) Xint[n-(i-1),] <- Xm[i,]

 intensity <- mean(X)
 x <- seq(0, max(X), 1)
 pexp <- dpois(x,lambda=intensity)
 XX <- table(X)
 pobs <- XX/sum(XX)

 X2 <- sum(((pobs-pexp)^2)/ pexp); X2
 df <- length(XX) - 2; df
 p.value <- 1 - pchisq(X2, df)

 panel4(size=7); par(pty="s")
 plot(pppset$x,pppset$y, xlab="x",ylab="y")
  #draw grid lines 
  for(i in 1:(n+1))  abline(h=ybrks[i],lty=2)
  for(i in 1:(n+1))  abline(v=xbrks[i],lty=2)
  title("Point Pattern",cex.main=0.7)
 # use transpose of Xm. We could have applied the im(.) command.  
 image(xbrks, ybrks, t(Xm), col=gray((10:0)/10),xlab="x",ylab="y")
  for(i in 1:(n+1))  abline(h=ybrks[i])
  for(i in 1:(n+1))  abline(v=xbrks[i])
  title("Intensity",cex.main=0.7)
  barplot(cbind(pobs,pexp),beside=T,ylab="Proportion"); 
  title("Counts Observed vs. Expected",cex.main=0.7)
 plot(ecdf(X),main="",ylab="F(Count)",xlab="Count",xlim=c(1,max(X))); points(x,ppois(x,lambda=intensity),col=1)
 title("ECDF",cex.main=0.7)


return(list(pppset=pppset, num.cells=n^2, Xint=Xint, chisq=X2, df=df, p.value=p.value, intensity=intensity))
}

# -------------------------
makeppp <- function(dataset){
# Reads a data frame and converts to ppp

#determine range
xr <- c( floor(min(dataset[,1])), ceiling(max(dataset[,1])) )
yr <- c( floor(min(dataset[,2])), ceiling(max(dataset[,2])) )

# make ppp
pppset <- ppp(dataset[,1],dataset[,2],xrange=xr,yrange=yr)

# describe it
summary(pppset)
return(pppset)
}

# -------------------------
nnGK.ppp <- function(dataset){

# make ppp
pppset <- makeppp(dataset)

# find n points
npoints <- length(dataset[,1])

G.pppset <- Gest(pppset, correction=c("none","km"))
K.pppset <- Kest(pppset, correction=c("none","iso"))

rG <- G.pppset$r[min(which(G.pppset$rs==1))+2]
rK <- max(K.pppset$r)

split.screen(c(2,1))

screen(1)
r <- G.pppset$r 
par(mar=c(4,4,1,.5),xaxs="r", yaxs="r")
  plot(G.pppset$raw~r, type="p", xlim=c(0,r[min(which(G.pppset$raw==1))]),xlab="Distance",ylab="Probability",cex=0.6)
  lines(G.pppset$theo~r, lty=1)
  lines(G.pppset$km~r, lty=2)
  legend("bottomright", c("Raw Uncorrected", "Theoretical Poisson","K-M Corrected"),
         pch=c(1,-1,-1), lty=c(-1,1,2), merge=T,cex=0.7)

split.screen(c(1,2),screen=2)
r <- K.pppset$r
screen(3)
  plot(K.pppset$un~r, type="p", xlim=c(0,max(r)),xlab="Distance",ylab="K(d)",cex=0.6)
  lines(K.pppset$theo ~r, lty=1)
  lines(K.pppset$iso~r, lty=2)
  legend("topleft", c("Raw Uncorrected", "Theoretical Poisson","Iso Corrected"),
         pch=c(1,-1,-1), lty=c(-1,1,2), merge=T,cex=0.7)
screen(4)
  plot(sqrt(K.pppset$un/pi)~r,  type="p", xlim=c(0,max(r)),ylim=c(0,max(r)),xlab="Distance",ylab="L(d)",cex=0.6)
  lines(sqrt(K.pppset$theo/pi) ~r,  lty=1)
  lines(sqrt(K.pppset$iso/pi)~r, lty=2)
  legend("topleft", c("Raw Uncorrected", "Theoretical Poisson","Iso Corrected"),
         pch=c(1,-1,-1), lty=c(-1,1,2), merge=T,cex=0.7)
}

# -------------------------
GKhat.env <- function(n, s, hat, stat, win){
 #
 dist <- unique(hat$r)

 # create place holder for s runs (rows) and distances (columns)
 hold <- matrix(0, s, length(dist))

 # make s Monte Carlo runs
 for(i in 1:s) {
 if(stat=="G") {hold[i, ] <- Gest(rpoispp(n,win=win),correction="km")$km
                xlabel <- "Mean Ghat";  ylabel <- "Ghat Empirical & Envelope"
		    leg.label <- "K-M"; xt <- hat$km; lim <-  1
               }  
 else          {hold[i, ] <- Kest(runifpoint(n,win=win))$iso
		    xlabel <- "Mean Khat";  ylabel <- "Khat Empirical & Envelope" 
		    leg.label <- "Iso"; xt <- hat$iso ; lim <- max(dist)
               }
 }

  #Calculate mean, max and min of all runs at each distance
  mn <- apply(hold,2, mean)
  Up <- apply(hold,2 , max)
  Down <- apply(hold,2, min)
  # combine by columns to plot
  y <- cbind(xt, mn, Down, Up) 

  # plot
  matplot(mn, y, type="l", lty=c(1,2,3,3), col=1, lwd=2,
          xlab=xlabel, ylab=ylabel, xlim=c(0,lim), c(0,lim))
  legend("topleft", legend=c(leg.label, "Mean", "Low&High"), lwd=2,lty=1:3, col=1)

 }

# -------------------------
nnGKenv.ppp <- function(dataset,nsim){

# make ppp
pppset <- makeppp(dataset)

# find n points
npoints <- length(dataset[,1])

G.pppset <- Gest(pppset)
K.pppset <- Kest(pppset)

panel2(size=7)
pppsetG.env <- GKhat.env(n=npoints, s=nsim, G.pppset, stat="G", win=pppset$window)
pppsetK.env <- GKhat.env(n=npoints, s=nsim, K.pppset, stat="K", win=pppset$window)

}

# -------------------------
vario <- function(dataset, num.lags, type='isotropic', theta, dtheta, maxdist){
 x <- names(dataset)[1]
 y <- names(dataset)[2] 
 v <- names(dataset)[3] 

 point.pair <- point(dataset)
 data.pair  <- pair(point.pair,num.lags=num.lags, type=type, theta, dtheta, maxdist=maxdist)
 dataset.v     <- est.variogram(point.pair,data.pair,v) 

panel4(size=7); par(cex.main=0.7)
plot.point(point.pair,v=v,legend.pos=2,pch=c(21:24),cex=0.7)
title("Dataset")
mtext(side=1,line=2,"x") 
mtext(side=2,line=2,"y")
plot(dataset.v)
spacecloud(point.pair,data.pair,v)
spacebox(point.pair,data.pair,v)
title("Boxplot by bin")
return(dataset.v)
}

# ------------------------
model.semivar.cov <- function(var, nlags, n0, c0, a){
# parameters
# var is empirical variogram, hmax is max of lag
# n0, c0, a: nugget, sill and range

# continuous lag
h <- seq(0,var$bins[nlags],var$bins[1]/100)

# Spherical model
g <- n0 + (c0-n0) * ( 3*h/(2*a) - h^3/(2*a^3) )
for(i in 1:length(g)) if(h[i] >=a) g[i] <- c0

# calculate covariance as sill minus semi-variance
c <- c0 - g
c[1] <- c0
g[1] <- 0
par(mfrow=c(2,1))
plot(h,g, type="n", ylim=c(0,max(var$classic/2)),xlab="Lag Distance h", ylab="Semi-variance gamma(h)")

# add empirical points note that we divide in half
points(var$bins,var$classic/2)

# plot the model
lines(h,g)

# plot the covariance
plot(h,c, type="l", ylim=c(0,1.1*c0),xlab="Lag Distance h", ylab="Covariance c(h)")
return()
}

# -----------------------------
# Function to put a matrix in a map for an image

img.map <- function(map){
 	image <- map
        nrows <- dim(map)[1]  
 	for(i in 1:nrows)
  	    image[nrows-(i-1),] <- map[i,]
        img <- t(image)
 return(img)
} 

