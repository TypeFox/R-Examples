
library(spatstat)
library(bayesm)

D <- data(longleaf)

x <- log((longleaf$x+.5)/201) - log(1-(longleaf$x+.5)/201)
y <- log((longleaf$y+.5)/201) - log(1-(longleaf$y+.5)/201)
m <- log(longleaf$marks-1.95)
Z <- cbind(x,y,m)
n <- nrow(Z)

lambda <- c(0,0,1)
kappa <- 0.01
nu <- 4
gam0 <- 4
psi0 <- diag(.1,3)
params <- c(lambda, kappa, nu, gam0, psi0)
N <- 200
mxd <- mix(Z, alpha=2, g0params=params,  N=N, niter=5, read=0, print=1)

nn <- 100
xx <- seq(1,199,length=(nx<-40))
mm <- seq(2.05,80, length=nn)
rr <- (rbind(c(100,100),c(150,150),c(150,50),c(25,150))+.5)/201
zm <- log(rr) - log(1-rr)
ZZ <- cbind(matrix(c(rep(zm[1,],nn), rep(zm[2,],nn), rep(zm[3,],nn), rep(zm[4,],nn)), ncol=2,byrow=TRUE), rep(log(mm-1.95),4))


dens <- function(prt){
  require(mvtnorm)

  # draw the relative weightings and expand with respect to G0
  probs <- rdirichlet(prt$n) 
  L <- 50
  v <- c(rbeta(L-1, 1, prt$n[1]), 1)
  for(l in 1:(L-1)){ v[l+1] <- v[l+1]*(1-sum(v[1:l])) }
  print(v[L])
  probs <- c(probs[1]*v, probs[-1])

  Omega <- matrix(as.numeric(prt[1,grep("S.",names(prt),fixed=TRUE)]), ncol=3)
  m <- (nrow(prt)-1)
  mu <- matrix(nrow=L+m, ncol=3)
  sigma <- array(dim=c(L+m,3,3))

  for(j in 1:L){
        sigma[j,,] <- solve(Bmix:::rwsh(nu, Omega))
        mu[j,] <- rmvnorm(1, lambda, sigma[j,,]/(kappa))
      }
  
  for(j in 1:m){    
    ### First draw all of the sigmas, etc...
    S <- matrix(as.numeric(prt[j+1,grep("S.",names(prt),fixed=TRUE)]), ncol=3)
    S[abs(S)<10^{-10}] <- 0
    zbar <- as.numeric(prt[j+1,grep("mean.",names(prt),fixed=TRUE)])
    nj <- prt$n[j+1]
    D <-  S + sum((lambda-zbar)^2)*kappa*nj/( kappa + nj )
    a <- as.numeric(prt[j+1,grep("a.",names(prt),fixed=TRUE)])
    sigma[j+L,,] <- solve(Bmix:::rwsh(nu+nj, Omega+D))
    mu[j+L,] <- rmvnorm(1, a , sigma[j+L,,]/(kappa+nj))
  }

  ## Build the densities
  zpdf <- matrix(rep(0,n), ncol=4)
  zzpdf <- rep(0,nn*4)
  zmarg <- rep(0,n)
  zzmarg <- rep(0,4)
  zcdf <- rep(0,n)
  zzcdf <- rep(0,nn*4)

  for(j in 1:(L+m)){
    for(i in 1:n){
      zmarg[i]= zmarg[i] + probs[j]*dmvnorm(Z[i,1:2], mu[j,1:2], sigma[j,1:2,1:2])
      zpdf[i]= zpdf[i] + probs[j]*dmvnorm(Z[i,], mu[j,], sigma[j,,])
    }
    for(i in 1:4){
      zzmarg[i]= zzmarg[i] + probs[j]*dmvnorm(ZZ[i*nn-1,1:2], mu[j,1:2], sigma[j,1:2,1:2])
    }
    for(i in 1:(4*nn)){
      zzpdf[i]= zzpdf[i] + probs[j]*dmvnorm(ZZ[i,], mu[j,], sigma[j,,])
    }
    print(j)
  }
  for(i in 1:n){ zcdf[i] <- zpdf[i]/zmarg[i] }
  for(k in 0:3){
    for(i in 1:nn){
      zzcdf[k*nn + i] <- zzpdf[k*nn + i]/zzmarg[k+1]
    }
  }
  return(list(zcdf=zcdf, zzcdf=zzcdf))
}


prts <- vector(mode="list", length=N)
for(i in 1:N) prts[[i]] <- particle(i, t=1, mxd)

post <- lapply(prts,dens)
zcdfs <- array( unlist(sapply(post,"[",1)), dim=c(n,n,N) )
zzcdfs <- array( unlist(sapply(post,"[",2)), dim=c(nn,4,N) )

pdf("pinesQQ.pdf",  paper="special", width=3, height=4.5, pointsize=8,family="Times")
par(mfrow=c(2,1), mai=c(.55,.55,.2,.1))
plot(longleaf$x, longleaf$y, cex=longleaf$marks/25,
     xlab="X1", ylab="X2", font.lab=3)
plot(seq(.01,1,length=100),seq(.01,1,length=100),type="l", col=grey(.8),
     lwd=2, xlab="Uniform Quantile", ylab="Estimated Quantile", cex=.5,   font.lab=3, cex.lab=1.2)
lines(1:n/(n+1), sort(apply(zcdfs, 1, mean)))
lines(1:n/(n+1), sort(apply(zcdfs, 1, quantile, .95)), lty=2)
lines(1:n/(n+1), sort(apply(zcdfs, 1, quantile, .05)), lty=2)
dev.off()


pdf("pineCND.pdf",  paper="special", width=4, height=3.2, pointsize=8,family="Times")
par(mfrow=c(2,2), mai=c(0,0,0,0), omi=c(.5,.5,.1,.1))
plot(mm, apply(zzcdfs[,1,],1,mean)/(mm-1.95), lty=1,type="l", yaxt="n", xaxt="n", ylab="n", xlab="n", ylim=c(0,0.1))
lines(mm, apply(zzcdfs[,1,],1,quantile,.95)/(mm-1.95), lty=2)
lines(mm, apply(zzcdfs[,1,],1,quantile,.05)/(mm-1.95), lty=2)
axis(side=2,  at=c(0,.05,.1,.15))
text(x=60, y=.08, labels="X = [100, 100]")
plot(mm, apply(zzcdfs[,2,],1,mean)/(mm-1.95), lty=1,type="l", yaxt="n", xaxt="n", ylab="n", xlab="n", ylim=c(0,0.1))
lines(mm, apply(zzcdfs[,2,],1,quantile,.95)/(mm-1.95), lty=2)
lines(mm, apply(zzcdfs[,2,],1,quantile,.05)/(mm-1.95), lty=2)
text(x=60, y=.08, labels="X = [150, 150]")
plot(mm, apply(zzcdfs[,3,],1,mean)/(mm-1.95), lty=1,type="l", yaxt="n", xaxt="n", ylab="n", xlab="n", ylim=c(0,0.1))
lines(mm, apply(zzcdfs[,3,],1,quantile,.95)/(mm-1.95), lty=2)
lines(mm, apply(zzcdfs[,3,],1,quantile,.05)/(mm-1.95), lty=2)
axis(side=2, at=c(0,.05,.1,.15))
axis(side=1, at=c(2,20,40,60,80))
text(x=60, y=.08, labels="X = [150, 50]")
plot(mm, apply(zzcdfs[,4,],1,mean)/(mm-1.95), lty=1,type="l", yaxt="n", xaxt="n", ylab="n", xlab="n", ylim=c(0,0.1))
lines(mm, apply(zzcdfs[,4,],1,quantile,.95)/(mm-1.95), lty=2)
lines(mm, apply(zzcdfs[,4,],1,quantile,.05)/(mm-1.95), lty=2)
axis(side=1, at=c(2,20,40,60,80))
text(x=60, y=.08, labels="X = [25, 150]")
mtext("Conditional Density", side=2, font=3, outer=TRUE, cex=1.1,line=3)
mtext("Diameter", side=1, font=3, outer=TRUE, cex=1.1,line=3)
dev.off()
