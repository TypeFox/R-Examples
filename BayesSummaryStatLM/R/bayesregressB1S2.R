## This method uses flat prior for beta, 1/sigmasq prior for sigmasq

bayesregressB1S2 <- function(xtx,xtx.inv,xty,yty, Tsamp.out,
                             sigmasq.init,numsamp.data)
{      
# define vectors and matrices

ytx<-t(xty)

n <- numsamp.data

Tsamp.out <- Tsamp.out+1

# num.predictors <- dim(xtx)[1]

betahat <- matrix(NA, nrow=Tsamp.out, ncol = dim(xtx)[1])

sigmasqhat <- rep(NA,Tsamp.out)

# set starting value for sigmasqhat

sigmasqhat[1] <- sigmasq.init

# posterior mean of betahat

betahat.mean <- xtx.inv%*%(xty)

# posterior variance of betahat

for (i in 2:Tsamp.out){

betahat[i,] <- rmvn(n=1, mu=betahat.mean, sigmasqhat[i-1]*xtx.inv)

# simulate sigmasqhat

sigmasqscale.pre <- (yty - t(betahat[i,]) %*% xty - 
ytx %*% betahat[i,] + t(betahat[i,]) %*% xtx %*% betahat[i,])

sigmasqhat[i] <- 1/rgamma(1,shape=n/2,scale=(sigmasqscale.pre/2)^(-1))

}  # end i

# remove starting value for sigmasqhat and NA for starting value of betahat

betahat    <- betahat[-1,]
sigmasqhat <- sigmasqhat[-1]

return(list("beta"=betahat,"sigmasq"=sigmasqhat))

}
