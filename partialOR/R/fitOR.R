# 12/31/2012 12:49:50 PM

fitOR <- function(dd) {
# ML-estimation of the F- and H-models
# dd: data frame,
# dd[,1:2] have names x, y; x,y are 0/1
# dd[, 4:(m+3)] are m numerical covariates

llikh <- function(ab,d2) {
# log-likelihood of the H model
# ab are model parameters (see the calling program partialORex.r)
  dz <- d2[,-1,drop=FALSE] 
  m <- dim(dz)[2]
  n <- dim(dz)[1]
  a <- ab[1:3]
  b <- matrix( ab[-(1:3)],nrow=2,byrow=TRUE)
  b <- rbind(b,b[1,]+b[2,])  # the homog.constraint
  ll <- 0
  for (i in 1:n) {
    h <- rep(0,4)
    for (j in 1:3) h[j+1] <- a[j] + drop( b[j,] %*% t(dz[i,]))
    hh <- sum(exp(h))
    ll <- ll + h[d2$xy[i]] - log(hh)
  }
  ll
} # end of llikh

dllikh <- function(ab,d2) {  
# 1st derivative of log-likelihood of the H model 
  dz <- d2[,-1,drop=FALSE]
  m <- dim(dz)[2]  
  n <- dim(dz)[1]
  a <- ab[1:3]
  b <- matrix(ab[-(1:3)],nrow=2,byrow=TRUE)
  b <- rbind(b,b[1,]+b[2,])  # b[,3] satisfies the homog.constraint
  dl <- rep(0,3+2*m)
  for (i in 1:n) {
    h <- rep(0,4)
    for (j in 1:3) h[j+1] <- a[j] + drop(b[j,] %*% t(dz[i,]))
    hh <- sum(exp(h))
    pp <- exp(h)/hh
    xy <- d2$xy[i]   
    for (j in 1:3) dl[j] <- dl[j] + (xy-1==j)*1 - pp[j+1] # dLL/daj
    for (k in 1:m) {
      k1 <- m + k
      ix <- ifelse ((xy==2 | xy==4),1,0)  # !
      iy <- ifelse ((xy==3 | xy==4),1,0)  # !
      dl[3+k]  <- dl[3+k] +(ix - pp[4]-pp[2])*dz[i,k]   # dLL/db1k
      dl[3+k1] <- dl[3+k1]+(iy - pp[4]-pp[3])*dz[i,k]   # dLL/db2k
    }
  }
  dl
} # end of dllikh

m <- dim(dd)[2]-2
covnam <- names(dd)[-(1:2)]
xy <- 1+2*dd$x+dd$y  # thus 1:4
ddxy <- cbind(xy,dd)
mmu  <- paste("xy ~",paste(covnam,collapse="+"),collapse="")
d2 <- ddxy[,c(-2,-3)]   # only the vars xy, z1,....,zm

# fit the null model:
fit0 <- multinom(xy~1, data=ddxy)   
# fit the full multinomial model :
sF <- summary(fitF <- multinom(formula(mmu), data=ddxy, Hess=TRUE))

# fit the H model using the F model parameters as start values :
start <- coef(sF)
start <- c(start[1:3,1],start[1,2:(m+1)],start[2,2:(m+1)])  
fitH <- optim(start, llikh, dllikh, d2, hessian=TRUE, method="BFGS", 
              control=list(fnscale=-1,maxit=200,trace=0) )  

cH <- matrix(NA,nrow=3,ncol=m+1)
cH[,1] <- fitH$par[1:3]              
for (j in 1:m) cH[1:2,1+j] <- fitH$par[-c(1:3)][(2*j-1):(2*j)]
cH[3,2:(m+1)] <- cH[1,2:(m+1)]+cH[2,2:(m+1)]
attributes(cH) <- attributes(fitF$coeff)              
fitH$coeff <- cH

return(list(fitH=fitH,fitF=sF,fit0=fit0))

} # end of fitOR

