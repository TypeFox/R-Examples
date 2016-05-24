qscatter_array<-function(variables,with.variables,data=NULL,x.lab="",y.lab="",
		main="Correlation Array",common.scales=TRUE,alpha=.25){
	arguments <- as.list(match.call()[-1])
	variables<-eval(substitute(variables),data,parent.frame())
	if(length(dim(variables))<1.5){
		variables<-d(variables)
		fn <- arguments$variables
		names(variables)<-if(is.call(fn)) format(fn) else as.character(fn)
	}
	if(missing(with.variables))
		with.variables <-variables
	else{
		with.variables<-eval(substitute(with.variables),data,parent.frame())
		if(length(dim(with.variables))<1.5){
			with.variables<-d(with.variables)
			fn <- arguments$with.variables
			names(with.variables)<-if(is.call(fn)) format(fn) else as.character(fn)
		}		
	}
	data1<-variables
	data2<-with.variables
	n<-nrow(data1)
	tmp<-rep(NA,n*ncol(data1)*ncol(data2))
	tmp.data<-data.frame(x=tmp,y=tmp,x.var=tmp,y.var=tmp)
	filled<-FALSE
	st<-1
	x<-1
	y<-1
	while(!filled){
		tmp.data[st:(st+n-1),1]<-data1[,x]
		tmp.data[st:(st+n-1),2]<-data2[,y]
		tmp.data[st:(st+n-1),3]<-names(data1)[x]
		tmp.data[st:(st+n-1),4]<-names(data2)[y]
		if(x>=ncol(data1) && y>=ncol(data2))
			filled<-TRUE
		if(y>=ncol(data2)){
			x<-x+1
			y<-1
		}else
			y<-y+1
		st<-st+n
	}
	p<-ggplot(tmp.data,aes(y,x))+geom_point(alpha=alpha)+facet_grid(x.var~y.var,scales="free")
	if(!common.scales)
		p<-p+facet_wrap(x.var~y.var,ncol=ncol(data2),scales="free")
	if(ncol(data1)==1 && y.lab=="")
		p<-p+ylab(names(data1)[1])
	else
		p<-p+ylab(y.lab)
	if(ncol(data2)==1 && x.lab=="")
		p<-p+xlab(names(data2)[1])
	else
		p<-p+xlab(x.lab)
	p+ggtitle(main)
}