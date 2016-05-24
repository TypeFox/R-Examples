plot.pd2bart = function(
   x,
   plquants =c(.05,.95), contour.color='white',
   justmedian=TRUE,
   ...
)
{
   pdquants = apply(x$fd,2,quantile,probs=c(plquants[1],.5,plquants[2]))
   qq <- vector('list',3)
   for (i in 1:3) 
      qq[[i]]  <- matrix(pdquants[i,],nrow=length(x$levs[[1]]))
   if(justmedian) {
     zlim = range(qq[[2]])
     vind = c(2)
   } else {
      par(mfrow=c(1,3))
      zlim = range(qq)
      vind = 1:3
   }
   for (i in vind) {
     image(x=x$levs[[1]],y=x$levs[[2]],qq[[i]],zlim=zlim,
       xlab=x$xlbs[1],ylab=x$xlbs[2],...)
     contour(x=x$levs[[1]],y=x$levs[[2]],qq[[i]],zlim=zlim,
       ,add=TRUE,method='edge',col=contour.color)
     title(main=c('Lower quantile','Median','Upper quantile')[i])
   }
}
