plotbranchmap<-function(bm,phi=55,theta=30)
{

persp(x=bm$level,y=bm$h,z=bm$z, 
xlab="level",ylab="h",zlab="excess mass",
ticktype="detailed", border=NA, shade=0.75,
col=bm$col,phi=phi,theta=theta) 

}


