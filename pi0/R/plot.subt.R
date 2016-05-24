plot.subt=function(x,y,rgl=TRUE,...)
{
    if(missing(y)) y=rgl else rgl=y
    n1=x[,'n1']
    n2=x[,'n2']
    y=x[,'f1']

    if(rgl) rgl=requireNamespace("rgl",quitely=TRUE)

    if(rgl){
        plot3d(n1,n2,y,col=4,size=3,xlab='n1',ylab='n2',zlab='pi0',...=...)
    }else{
        warning("package rgl not found/used. scatterplot3d is used.")
        #library("scatterplot3d")
        scatterplot3d(n1,n2,y,color=4,pch=19,xlab='n1',ylab='n2',zlab='pi0',
            xlim=c(1,max(n1)),ylim=c(1,max(n2)),zlim=range(y),...=...)
    }
    invisible(NULL)
}
