gng.plot.fit <-
function (data, obj, resolution=100,breaks=100,legpos=NULL,xlim=NULL,
  main=NULL,...)
{  
   obs <- unlist(data);
   if(is.null(xlim)){ xlim=range(obs);}
   if(is.null(main)){main="Goodness of Fit";}             
   hist(obs,freq=FALSE,breaks=breaks,xlim=xlim,main=main,...);
   gng.plot.mix(obj,resolution=resolution,col='black',lwd=3,new.plot=FALSE,...);
   gng.plot.comp(data,obj,new.plot=FALSE,xlim=xlim,legpos=legpos,lwd=2,...);
}

