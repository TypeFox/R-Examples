########## R-function: plotControl ##########

# For controlling the plot parameters for plot.spm()

# Last changed: 18 JAN 2005

plotControl <- function(plot.it=NULL,drv=0,
                         se=TRUE,shade=TRUE,
                         bty=NULL,main=NULL,xlab=NULL,
                         ylab=NULL,xlim=NULL,ylim=NULL,grid.size=NULL,
                         lty=NULL,lwd=NULL,col=NULL,
                         se.lty=NULL,se.lwd=NULL,se.col=NULL,
                         shade.col=NULL,rug.col=NULL,
                         jitter.rug=FALSE,zero.line=NULL,
                         plot.image=TRUE,
                         image.col=blend.col("cyan","red"),
                         image.bg=NULL,image.bty=NULL,
                         image.main=NULL,image.xlab=NULL,
                         image.ylab=NULL,image.zlab="",
                         image.xlim=NULL,image.ylim=NULL,
                         image.zlim=NULL,
                         image.grid.size=NULL,bdry=NULL,
                         add.legend=TRUE,leg.loc=NULL,
                         leg.dim=NULL,image.zlab.col=NULL)
                                                         
   return(list(plot.it=plot.it,drv=drv,se=se,shade=shade,
               bty=bty,main=main,xlab=xlab,ylab=ylab,
               xlim=xlim,ylim=ylim,grid.size=grid.size,
               lty=lty,lwd=lwd,col=col,
               se.lty=se.lty,se.lwd=se.lwd,se.col=se.col,
               shade.col=shade.col,rug.col=rug.col,
               jitter.rug=jitter.rug,
               zero.line=zero.line,plot.image=plot.image,
               image.col=image.col,image.bg=image.bg,
               image.bty=image.bty,image.main=image.main,
               image.xlab=image.xlab,image.ylab=image.ylab,
               image.zlab=image.zlab,image.xlim=image.xlim,
               image.ylim=image.ylim,image.zlim=image.zlim,
               image.grid.size=image.grid.size,bdry=bdry,
               add.legend=add.legend,
               leg.loc=leg.loc,leg.dim=leg.dim,
               image.zlab.col=image.zlab.col))
                         
######### End of plotControl ##########












