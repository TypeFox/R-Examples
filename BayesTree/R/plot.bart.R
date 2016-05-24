plot.bart = function(
   x,
   plquants=c(.05,.95), cols =c('blue','black'),
   ...
)
{
   if('sigma' %in% names(x)) {
   par(mfrow=c(1,2))
   plot(c(x$first.sigma,x$sigma),col=rep(c('red','black'),
     c(length(x$first.sigma),length(x$sigma))),ylab='sigma',...)
   } else {
     par(mfrow=c(1,1))
   }

  if('sigma' %in% names(x)) {
     ql <- apply(x$yhat.train,2,quantile,probs=plquants[1])
     qm <- apply(x$yhat.train,2,quantile,probs=.5)
     qu <- apply(x$yhat.train,2,quantile,probs=plquants[2])
     plot(x$y,qm,ylim=range(ql,qu),xlab='y',ylab=
      'posterior interval for E(Y|x)',...)
     for (i in 1:length(qm))
       lines(rep(x$y[i],2),c(ql[i],qu[i]),col=cols[1])
     abline(0,1,lty=2,col=cols[2])
  } else {
     pdrs = pnorm(x$yhat.train) #draws of p(Y=1 | x)
     ql <- apply(pdrs,2,quantile,probs=plquants[1])
     qm <- apply(pdrs,2,quantile,probs=.5)
     qu <- apply(pdrs,2,quantile,probs=plquants[2])
     plot(qm,qm,ylim=range(ql,qu),xlab='meadian of p',ylab=
      'posterior interval for P(Y=1|x)',...)
     for (i in 1:length(qm))
       lines(rep(qm[i],2),c(ql[i],qu[i]),col=cols[1])
     abline(0,1,lty=2,col=cols[2])
  }
}
