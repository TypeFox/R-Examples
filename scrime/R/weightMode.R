`weightMode` <-
function(w,v,num=TRUE,samp=FALSE){
	tab<-by(w,v,sum,na.rm=TRUE)
	if(!samp)
		out<-names(tab)[tab==max(tab)]
	else
		out<-if(length(tab)==1) names(tab) else sample(names(tab),1,prob=tab)
	if(length(out)>1)
		out<-sample(out,1)
	if(num)
		out<-as.numeric(out)
	out
}

