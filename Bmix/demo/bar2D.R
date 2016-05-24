###############
# 2D dynamic density estimation example
###############

# Simulate the data

library(mvtnorm)
nz <- 20
z <- array(dim=c(nz,3,10))
y <- matrix(ncol=2, nrow=nz)
Z <- c()
for(i in 1:10)
{  x <- rbinom(nz,1, prob=(i+3)/14)
   if(sum(1-x)==0){ x[nz] <- 0 }  #just to stop errors
   if(sum(x)==0){ x[nz] <- 1 }
   y<- rbind(rmvnorm(sum(1-x),c(0,0),diag(.5,2)),
            rmvnorm(sum(x),c(5,0),diag(2,2))) 
   z[,,i] <- cbind(sort(x),y)[sample(1:nz, nz),]
   Z<- rbind(Z,z[,2:3,i]) }

N <- 100; alpha <- 4; rho <- 0.8
params <- c(c(2,0), #gamma
            0.1, #kappa
            3, #nu
            3, #gam0
            diag(.5,2) #psi0
            )
times <- c()
for(i in 0:9){ times <- c(times, rep(i, nz)) }

# Filter BAR stick-breaking mixture with rho=.8
mxdyn <- mix(Z, alpha=alpha, g0params=params, times=times, rho=rho,
	 		 	cat=0, N=N, niter=0, read=0, print=1) 

# Read densities in from text
prts <- vector(mode="list", length=0)
for(t in 1:10){
  prt <-  vector(mode="list", length=N) 
  for(i in 1:N) prt[[i]] <- particle(i, mxdyn, t, rho)
  prts <- cbind(prts, prt) }

# Grid for prediction
xx <- seq(-2,10,length=(nx<-40))
yy <- seq(-5,5,length=(ny<-40))
zz <- expand.grid(xx,yy)

# function to extract various densities from each particle
dens <- function(prt)
  { require(mvtnorm)
    pdf <- matrix(rep(0,nx*ny), ncol=ny,nrow=nx)
    for(j in 1:nrow(prt))
      { pdf <- pdf + prt$p[j]*dmvt(t(t(zz)-as.numeric(prt[j,grep("a.",names(prt),fixed=TRUE)])),
                                   sigma=matrix(as.numeric(prt[j,grep("B.",names(prt),fixed=TRUE)]),
                                     ncol=2), df = prt$c[j], log=FALSE) }
    return(pdf) }

# Extract mean pdfs
post <- lapply(prts,dens)
pdfs <- array( unlist(post), dim=c(nx,ny,N,10) )

# Plot Filtered Densities
cols <- c(5,6)
par(mfrow=c(2,5), mai=c(0,0,0,0), omi=c(.2,.2,0.05,0.05) )
for(j in 0:9){
  plot(Z[j*nz+1:nz,], col=cols[z[,1,j+1]+1],
       pch=20, xaxt="n", yaxt="n", xlim=range(xx), ylim=range(yy))
  contour(xx,yy, apply(pdfs[,,,j+1], c(1,2), mean),
          add=TRUE, levels=c(1:20)*.01, drawlabels=FALSE)
  if(j==0 || j==5){ axis(2, cex.axis=1.3) }
  if(j>4){ axis(1, cex.axis=1.3, at=c(0,2,4,6,8)) }
}

# Run a simple DP mixture model for the time t=10 observations
mxd <- mix(Z[9*nz+1:nz,], alpha=alpha, g0params=params, cat=0, N=N, niter=10,  print=1)
prt <-  vector(mode="list", length=N) 
for(i in 1:N) prt[[i]] <- particle(i, mxd, 1)
pf <- array( unlist(lapply(prt, dens)), dim=c(nx,ny,N) )
rl <- readline("press RETURN to continue: ")

# Comparison between filtered and independent fit
par(mfrow=c(1,3), mai=c(0,0,.2,.1), omi=c(.2,.2,.1,.1))
contour(xx,yy, apply(pdfs[,,,10], c(1,2), mean), xlim=range(xx), ylim=range(yy),
          levels=c(1:20)*.01, drawlabels=FALSE, xaxt="n", yaxt="n", main="Filtered AR Fit", cex.main=1.2)
points(Z[181:200,], col=grey(.5), pch=20)
axis(1, cex.axis=1.3, at=c(0,2,4,6,8))
axis(2, cex.axis=1.3)
contour(xx,yy, ((1/14)*matrix(dmvnorm(zz, c(0,0),diag(.5,2)), ncol=ny)
            + (13/14)*matrix(dmvnorm(zz,c(5,0),diag(2,2)), ncol=ny)),
        xlim=range(xx), ylim=range(yy),
        levels=c(1:20)*.01, drawlabels=FALSE, xaxt="n", yaxt="n", main="The Truth", cex.main=1.2)
points(Z[181:200,], col=grey(.5), pch=20)
axis(1, cex.axis=1.3, at=c(0,2,4,6,8))
contour(xx,yy, apply(pf, c(1,2), mean), xlim=range(xx), ylim=range(yy),
          levels=c(1:20)*.01, drawlabels=FALSE, xaxt="n", yaxt="n", main="Independent Fit", cex.main=1.2)
axis(1, cex.axis=1.3, at=c(0,2,4,6,8))
points(Z[181:200,], col=grey(.5), pch=20)

