`modeDist` <-
function(x,num=TRUE){
	tab<-table(x)
	out<-names(tab)[tab==max(tab)]
	if(length(out)>1)
		out<-sample(out,1)
	if(num)
		out<-as.numeric(out)
	out
}

