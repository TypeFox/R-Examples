plot.pdbart = function(
   x,
   xind = 1:length(x$fd),
   plquants =c(.05,.95),cols=c('black','blue'),
   ...
)
{
   rgy = range(x$fd)
   for(i in xind) {
         tsum = apply(x$fd[[i]],2,quantile,probs=c(plquants[1],.5,plquants[2]))
         plot(range(x$levs[[i]]),rgy,type='n',xlab=x$xlbs[i],ylab='partial-dependence',...)
         lines(x$levs[[i]],tsum[2,],col=cols[1],type='b')
         lines(x$levs[[i]],tsum[1,],col=cols[2],type='b')
         lines(x$levs[[i]],tsum[3,],col=cols[2],type='b')
   }
}
