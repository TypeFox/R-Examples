## This method uses normal prior for beta with unknown mean and variance, 
## Inverse gamma prior for sigmasq
## Rinv and Vinv are SPD matrices

bayesregressB3S1 <- function(xtx, xty, yty, numsamp.data,
                             eta       = rep(1.0, dim(xtx)[1]),
                             Rinv      = diag(1.0, dim(xtx)[1]),                             
                             mu.init   = rep(1.0, dim(xtx)[1]),
                             lambda    = dim(xtx)[1],
                             Vinv      = diag(1.0, dim(xtx)[1]),                             
                             Cinv.init = diag(1.0, dim(xtx)[1]),
                             inv.gamma.a = 1.0, 
                             inv.gamma.b = 1.0,
                             sigmasq.init=1,                             
                             Tsamp.out)
{
# define vectors and matrices

ytx<-t(xty)

n <- numsamp.data

Tsamp.out <- Tsamp.out+1

# num.predictors = dim(xtx)[1]

betahat <- matrix(NA,nrow=Tsamp.out,ncol=dim(xtx)[1])

sigmasqhat <- rep(NA,Tsamp.out)

muhat <- matrix(NA,nrow=Tsamp.out,ncol=dim(xtx)[1])

Cinvhat <- array(NA,c(Tsamp.out,dim(xtx)[1],dim(xtx)[1]))

# set starting value for sigmasqhat, muhat, Cinvhat

sigmasqhat[1] <- sigmasq.init

muhat[1,] <- mu.init

Cinvhat[1,,] <- Cinv.init

# posterior mean of betahat

for (i in 2:Tsamp.out){

betahat.var <- chol2inv(chol((Cinvhat[(i-1),,] + (1/sigmasqhat[i-1])*(xtx))))
betahat.mean.pre1 <- Cinvhat[(i-1),,] %*% muhat[(i-1),] +
 (1/sigmasqhat[i-1])* xty
betahat.mean <- betahat.var %*% betahat.mean.pre1

betahat[i,] <- rmvn(n=1,mu = betahat.mean,sigma = betahat.var)

muhat.var <- chol2inv(chol((Rinv+Cinvhat[(i-1),,])))
muhat.mean.pre1 <- (Cinvhat[(i-1),,] %*% betahat[i,] + Rinv %*% eta)
muhat.mean <- muhat.var %*% muhat.mean.pre1

muhat[i,] <- rmvn(n=1,mu = muhat.mean,sigma = muhat.var)

#Cinvhat[i,,] <- rWISHART(n=1,df=(lambda+1),chol2inv(chol((Vinv+(betahat[i,]-muhat[i,]) %*% t(betahat[i,]-muhat[i,])))))

Cinvhat[i,,] <- rWishart(n=1,df=(lambda+1),chol2inv(chol((Vinv+(betahat[i,]-muhat[i,]) %*% t(betahat[i,]-muhat[i,])))))

# simulate sigmasqhat

# sigmasqhat.pre

sigmasqscale.pre <- (yty - t(betahat[i,]) %*% xty - 
ytx %*% betahat[i,] + t(betahat[i,]) %*% xtx %*% betahat[i,])


sigmasqhat[i] <- 1/rgamma(1,shape=n/2+inv.gamma.a,scale=(sigmasqscale.pre/2+(1/inv.gamma.b))^(-1))

}  # end i

# remove initial values (some are NA's)

betahat    <- betahat[-1,]
sigmasqhat <- sigmasqhat[-1]
muhat      <- muhat[-1,]
Cinvhat    <- Cinvhat[-1,,]

return(list("beta"=betahat,"sigmasq"=sigmasqhat,"mu"=muhat,"Cinv"=Cinvhat))

}






