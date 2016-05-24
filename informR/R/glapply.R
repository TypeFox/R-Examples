`glapply` <-
function(x,id,FUN,regroup=TRUE,...){
	   x.l<-split(x,id)
	   fx.l<-lapply(x.l,FUN=FUN,...)
	   if(regroup){unsplit(fx.l,id)}
	   else fx.l
	}
