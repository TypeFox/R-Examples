plot.extrpi0=function(x,y,rgl=TRUE,...)
{
    if(missing(y)) y=rgl else rgl=y
    subt.data=attr(x,'subt.data')
    n1=subt.data[,'n1']
    n2=subt.data[,'n2']
    y=subt.data[,'f1']
    parm=attr(x,'par')
    f=function(fn1=n1,fn2=n2) parm[3]*((fn1+fn2)/(fn1+fn2+fn1*fn2*parm[2]))^parm[4]+parm[1]

    nn1=seq(1,2*max(n1),length=50)
    nn2=seq(1,2*max(n1),length=50)
    plane=outer(nn1,nn2,f)
    if(rgl) rgl=requireNamespace("rgl",quitely=TRUE)

    if(rgl){
        persp3d(nn1,nn2,plane,xlab='n1',ylab='n2',zlab='pi0',
                color='red',zlim=c(min(x[1],y),max(plane,y)),...)
        plot3d(n1,n2,y,add=TRUE,col=4,size=3)
    }else{
        warning("package rgl not found/used. Only static plot is generated")
        res=persp(nn1,nn2,plane,xlab='n1',ylab='n2',zlab='pi0',
                col='red',zlim=c(min(x[1],y),max(plane,y)),
                shade=.75,theta=50,phi=20,lphi=50,ltheta=-20,ticktype="detailed",...)
        points(trans3d(n1,n2,y,res),col=4,pch=19)
    }
    invisible(NULL)
}
