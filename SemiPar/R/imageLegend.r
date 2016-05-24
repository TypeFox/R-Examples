########## R function: imageLegend ##########

# For approximately mimicking the S-PLUS function image.legend()

# Last changed: 18 JAN 2005

imageLegend <- function(zlim,col,top.left.loc,dims,zlab,label.col)
{
   x.vec <- seq(top.left.loc[1],top.left.loc[1]+dims[1],length=1001)
   y.vec <- c(top.left.loc[2]-dims[2],top.left.loc[2])

   z.mat <- as.matrix(seq(zlim[1],zlim[2],length=1000))

   image(x.vec,y.vec,z.mat,col=col,add=T)

   lines(c(min(x.vec),max(x.vec)),rep(min(y.vec),2),col=label.col)
   lines(c(min(x.vec),max(x.vec)),rep(max(y.vec),2),col=label.col)
   lines(rep(min(x.vec),2),c(min(y.vec),max(y.vec)),col=label.col)
   lines(rep(max(x.vec),2),c(min(y.vec),max(y.vec)),col=label.col)

   # Add the variable name

   text(top.left.loc[1]+0.5*dims[1],top.left.loc[2]+0.35*dims[2],
        zlab,col=label.col)

   # Add the axis

   axs.ht <- y.vec[1]-0.3*dims[2]

   lines(c(min(x.vec),max(x.vec)),rep(axs.ht,2),col=label.col)

   tick.nums <- pretty(z.mat)

   tick.nums <- tick.nums[tick.nums<=max(z.mat)]
   tick.nums <- tick.nums[tick.nums>=min(z.mat)]

   zt.low <- top.left.loc[1]
   zt.upp <- top.left.loc[1]+dims[1]

   sca.fac <- (zt.upp-zt.low)/(zlim[2]-zlim[1])

   t.tick.nums <- zt.low + sca.fac*(tick.nums-zlim[1])

   for (i in 1:length(tick.nums))
   {
      lines(rep(t.tick.nums[i],2),c(axs.ht,axs.ht+0.2*dims[2]),col=label.col) 
      text(t.tick.nums[i],axs.ht-0.35*dims[2],
           as.character(tick.nums[i]),col=label.col)
   }
}

########## End of imageLegend ##########
