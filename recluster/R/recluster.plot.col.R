recluster.plot.col<-function (mat,cext=0.3,cex=1,cex.axis=0.7,cex.lab=0.8,pch=16,text=TRUE,add=F,xlim=NULL,ylim=NULL,ylab="Axis 2",xlab="Axis 1",...) 
{
        if (diff(range(mat[,1]))>diff(range(mat[,2]))){	  
                if (is.null(xlim)) {xlim=range(mat[,1])}
                if (is.null(ylim)) {
                      y=(diff(range(mat[,1]))-diff(range(mat[,2])))/2
                      ylim<-range(mat[,2])+c(-y,y)}
                      }
        if (!add){plot(mat[,1],mat[,2],col=rgb(mat[,3],mat[,4],mat[,5], maxColorValue = 255),cex=cex,xlim=xlim, ylim=ylim,pch=pch,ylab=ylab,xlab=xlab,...)}
        if (add){points(mat[,1],mat[,2],col=rgb(mat[,3],mat[,4],mat[,5], maxColorValue = 255),cex=cex,xlim=xlim, ylim=ylim,pch=pch,...)}
        if (text){text(mat[,1],mat[,2]-0.03,rownames(mat),cex=cext)}
}
