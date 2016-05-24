FM.chi.pvalue <-
function(x){
	y=log(x)
	y=-2*sum(y)
	 p.value=pchisq(y,df=2*length(x),lower.tail=FALSE)
	list("FMstatistic"=y, "FMpvalue"=p.value)
  }
