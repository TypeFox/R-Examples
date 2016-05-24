#####
# Example marginal likelihood calculations for alpha
######

library(MASS)

N <- 200 

alpha <- seq(0.1,15, length=16)


### Galaxy data
data(galaxies) 
gparams <- c(0, #gamma
            1/2, #kappa
            2, #nu
            2, #gam0
            1/2 #psi0
            )

x <- ((galaxies-mean(galaxies))/sd(galaxies))[sample(1:length(galaxies))]
gnum <- glhd <- c()
for(i in 1:length(alpha)){
  out <- mix(x, alpha=alpha[i], g0params=gparams, N=N)
  gnum <- c(gnum, out$m[length(x)])
  glhd <- c(glhd, out$logprob[length(x)])
}


### Animals data (residuals)
data(Animals)
attach(Animals)
aparams <- c(0,1/2,2,2,1/2)
reg <- summary(lin <- lm(log(brain) ~ log(body)))
r <- (reg$residuals/reg$sigma)[sample(1:nrow(Animals))]
anum <- alhd <- c()
for(i in 1:length(alpha)){
  out <- mix(r, alpha=alpha[i], g0params=aparams, N=N)
  anum <- c(anum, out$m[length(r)])
  alhd <- c(alhd,  out$logprob[length(r)])
}


### Motorcycle data
data(mcycle)
mm <- apply(mcycle,2,mean)
msd <- apply(mcycle,2,sd)
Z <- t( (t(mcycle)-mm)/msd )[sample(1:nrow(mcycle)),]
mparams <- c(c(0,0), 1/5, 3, 3, diag(1/5,2) )  #new params for 2D
mnum <- mlhd <- c()
for(i in 1:length(alpha)){
  out <- mix(Z, alpha=alpha[i], g0params=mparams, N=N)
  mnum <- c(mnum, out$m[nrow(Z)])
  mlhd <- c(mlhd,  out$logprob[nrow(Z)])
}

### Check that the models are doing what we think they are...


# Grids and predictive density functions

xx <- seq(-3,3,length=(nn<-40))
zz <- expand.grid(xx,xx)
gpdf <- rep(0,nn)
apdf <- rep(0,nn)
mpdf <- matrix(rep(0,nn^2), ncol=nn,nrow=nn)

dens1 <- function(prt){
  pdf <- rep(0,nn)
  for(j in 1:nrow(prt))
    { pdf <- pdf + prt$p[j]*dt( (xx-prt[j,]$a.1)/sqrt(prt[j,]$B.1),
                               df = prt$c[j])/sqrt(prt[j,]$B.1) }
  return(pdf) }

dens2 <- function(prt){
  require(mvtnorm)
  pdf <- matrix(rep(0,nn^2), ncol=nn,nrow=nn)
  for(j in 1:nrow(prt)){
    pdf <- pdf + prt$p[j]*dmvt(t(t(zz)-as.numeric(prt[j,grep("a.",names(prt),fixed=TRUE)])),
                               sigma=matrix(as.numeric(prt[j,grep("B.",names(prt),fixed=TRUE)]), ncol=2),
                               df = prt$c[j], log=FALSE) }
  return(pdf) }

# take the maximum likelihood alphas

galpha <- alpha[which.max(glhd)]
aalpha <- alpha[which.max(alhd)]
malpha <- alpha[which.max(mlhd)]


out <-  mix(x, alpha=galpha, g0params=gparams, N=N, print=TRUE)
for(i in 1:N){ gpdf = gpdf + dens1(particle(i, out, 1))/N }

out <-  mix(r, alpha=galpha, g0params=aparams, N=N, print=TRUE)
for(i in 1:N){ apdf = apdf + dens1(particle(i, out, 1))/N }

out <-  mix(Z, alpha=galpha, g0params=mparams, N=N, print=TRUE)
for(i in 1:N){ mpdf = mpdf + dens2(particle(i, out, 1))/N }


# plot it all
par(mai=c(.6,.6,0.3,0.1), mfrow=c(3,4))
hist(galaxies, col=8, main="galaxies data")
plot(alpha, gnum, type="l", lwd=1.5,
     xlab="alpha", ylab="mean number of components", main="galaxies M")
plot(alpha, glhd,type="l", lwd=1.5,
     xlab="alpha", ylab="log marginal likelihood", main="galaxies LHD")
plot(xx, gpdf, type="l", main=paste("alpha =", round(galpha,0)))
points(x, rep(0, length(x)), pch="*", col=8, cex=1.5)
plot(log(Animals), main="animals data", pch=20)
lines(log(body), lin$fitted.values)
plot(alpha, anum, type="l", lwd=1.5,
     xlab="alpha", ylab="mean number of components", main="Animals M")
plot(alpha, alhd, type="l", lwd=1.5,
     xlab="alpha", ylab="log marginal likelihood", main="Animals LHD")
plot(xx, apdf, type="l", main=paste("alpha =", round(aalpha,0)))
points(r, rep(0.01, length(r)), pch="*", col=8, cex=1.5)
plot(mcycle, main="motorcycle data", pch=20)
plot(alpha, mnum, type="l", lwd=1.5,
     xlab="alpha", ylab="mean number of components", main="mcycle M")
plot(alpha, mlhd, type="l", lwd=1.5,
     xlab="alpha", ylab="log marginal likelihood", main="mcycle LHD")
plot(Z, pch=20, col=8,  main=paste("alpha =", round(malpha,0))) 
contour(xx, xx, mpdf, add=TRUE, drawlabels=FALSE)
