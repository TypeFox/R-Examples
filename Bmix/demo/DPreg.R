######
# DP regression example, with both categorical and continuous covariate
######

# Simulate the data
x <- rnorm(n <- 500, 0, 1)
h <-  function(x){ 0.3 + 0.4*x + 0.5*sin(2.7*x) + 1.1*(1 + x^2 )^(-1) }
m <- rbinom(n,1,pnorm(x))
Z <- cbind(x, y <- h(x)+rnorm(n,0,c(.25,.5)[m+1]), m)

alpha <- 2
g0params <- c(c(0,0), #gamma
            0.01, #kappa
            3, #nu
            3, #gam0
            diag(.1,2) #psi0
            )

# Filter with 100 particles, and follow up with 10 Gibbs iterations
mxd <- mix(Z, alpha, g0params, cat=1, N=100, niter=10, print=1)
attach(mxd)

# Grid for plotting
xx <- seq(-3,3,length=(nx<-40))
yy <- seq(-1,3,length=(ny<-40))
zz <- expand.grid(xx,yy)

# function to extract various densities from each particle
dens <- function(prt){
  require(mvtnorm)
  pdf <- matrix(rep(0,nx*ny), ncol=ny,nrow=nx)
  mdf <- rep(0,nx)
  pf <- rep(0,nx)
  for(j in 1:nrow(prt)){
    pdf <- pdf + prt$p[j]*dmvt(t(t(zz)-as.numeric(prt[j,grep("a.",names(prt),fixed=TRUE)])),
                               sigma=matrix(as.numeric(prt[j,grep("B.",names(prt),fixed=TRUE)]), ncol=dim),
                               df = prt$c[j], log=FALSE)
    tmp <- prt$p[j]*dmvt(matrix(xx-prt$a.1[j]), sigma=matrix(prt$B.1[j]), df=prt$c[j], log=FALSE)
    mdf <- mdf + tmp
    if(j==1){ pf <- pf + tmp*prt$counts.1.2[1]/sum(prt[1,grep("counts", names(prt))]) }
    else{ pf <- pf + tmp*sum(prt$counts.1.2[c(1,j)])/sum(prt[c(1,j),grep("counts", names(prt))])}
  }
  # note that this is a quick but sloppy way to do conditioning;
  # see Taddy & Kottas 2010 JBES, as well as the pines.R demo
  cdf <- pdf/matrix(mdf,nrow=nx,ncol=ny)
  pf <- pf/mdf
  return(list(pdf=pdf,mdf=mdf,cdf=cdf,pf=pf))
}

# Read in particles from the text files
prts <- vector(mode="list", length=N)
for(i in 1:N) prts[[i]] <- particle(i, mxd, 1)

# Extract various densities
post <- lapply(prts,dens)
pdfs <- array( unlist(sapply(post,"[",1)), dim=c(nx,ny,N) )
cdfs <- array( unlist(sapply(post,"[",3)), dim=c(nx,ny,N) )
pfs <- matrix( unlist(sapply(post,"[",4)), nrow=nx, ncol=N )

# Plot 
par(mfrow=c(1,3))
plot(Z, pch=20,  col=8, main="joint pdf")
contour(xx,yy, apply(pdfs,c(1,2),mean),add=TRUE)
plot(Z, pch=20,  col=8, main="conditional pdf")
contour(xx,yy, apply(cdfs,c(1,2),mean),add=TRUE)
plot(xx, pnorm(xx),type="l", col=8, main="p(m=1)")
lines(xx, apply(pfs, 1, mean))
lines(xx, apply(pfs, 1, quantile, .95))
lines(xx, apply(pfs, 1, quantile, .05))
