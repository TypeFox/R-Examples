fracdev=function(U,y,alpha,beta,intercept,degrees,family){
  y=drop(y)# in case it was entered as a one-column matrix
  p=dim(alpha)[1]
  offset=c(1,1+cumsum(degrees)[-p])
  beta[offset,]=beta[offset,]+alpha
  fitmat=scale(U%*%beta,-intercept,FALSE)
  switch(family,
         gaussian={ nulldev=sum( (y-mean(y))^2)
                    dev=apply( (y-fitmat)^2,2,sum)
                    list(nulldev=nulldev,dev.ratio=1-dev/nulldev)
                  },
         binomial={ mu=mean(y);nulldev=-2*sum (y*log(mu)+(1-y)*log(1-mu))
                    pmat=1/(1+exp(-fitmat))
                    dev=-2*apply(y*log(pmat)+(1-y)*log(1-pmat),2, sum)
                    list(nulldev=nulldev,dev.ratio=1-dev/nulldev)
                  }
         )
}
