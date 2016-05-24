 dea_est<-function(xtab,ytab,x,type="dea")
 {
  stopifnot(type%in%c("dea","fdh","lfdh"),length(xtab)==length(ytab))
 
  if(type=="dea")
   {return((dea(x,rep(median(ytab),length(x)),XREF=xtab,YREF=ytab,RTS=1,ORIENTATION="out")$eff)*median(ytab))}
    else
    {if(type=="fdh")
     {return((dea(x,rep(median(ytab),length(x)),XREF=xtab,YREF=ytab,RTS=0,ORIENTATION="out")$eff)*median(ytab))
     }
     else
     {fdh_index<-which(dea(xtab,ytab,RTS=0,ORIENTATION="out")$eff==1)
      itp<-approx(x=xtab[fdh_index], y = ytab[fdh_index], x, method = "linear",yleft=min(ytab[fdh_index]),yright=max(ytab[fdh_index]))
      return(itp$y)
     }
   }
 }
 
 