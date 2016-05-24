## This method uses normal prior for beta, inverse gamma prior for sigmasq

bayesregressB2S1 <- function(xtx,xty,yty,numsamp.data,
                             beta.prior.mean = rep(0,dim(xtx)[1]),
                             beta.prior.var  = diag(1.0,dim(xtx)[1]),
                             beta.prior.var.inv = chol2inv(chol(beta.prior.var)),
                             inv.gamma.a = 1.0, inv.gamma.b = 1.0,
                             sigmasq.init = 1.0,
                             Tsamp.out)
{      
# define vectors and matrices

ytx<-t(xty)

n <- numsamp.data

Tsamp.out <- Tsamp.out+1

# num.predictors <- dim(xtx)[1]

betahat <- matrix( NA, nrow=Tsamp.out, ncol=dim(xtx)[1] )

sigmasqhat <- rep(NA,Tsamp.out)

# set starting value for sigmasqhat

sigmasqhat[1] <- sigmasq.init

betahat.pre1 <- beta.prior.var.inv
betahat.pre2 <- betahat.pre1 %*% beta.prior.mean

# posterior variance of betahat

for (i in 2:Tsamp.out){

# use cholesky for matrix inverse

betahat.var <- chol2inv(chol((betahat.pre1+(1/sigmasqhat[i-1]) * xtx)))

betahat.mean <- betahat.var %*% (betahat.pre2 + (1/sigmasqhat[i-1]) * xty)

betahat[i,] <- rmvn(n=1,mu = betahat.mean,sigma = betahat.var)

# simulate sigmasqhat

sigmasqscale.pre <- (yty - t(betahat[i,]) %*% xty - 
ytx %*% betahat[i,] + t(betahat[i,]) %*% xtx %*% betahat[i,])

sigmasqhat[i] <- 1/rgamma(1,shape=n/2+inv.gamma.a,scale=(sigmasqscale.pre/2+(1/inv.gamma.b))^(-1))

}  # end i

# remove starting value for sigmasqhat and NA for starting value of betahat

betahat    <- betahat[-1,]
sigmasqhat <- sigmasqhat[-1]

return(list("beta"=betahat,"sigmasq"=sigmasqhat))

}