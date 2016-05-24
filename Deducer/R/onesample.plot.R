onesample.plot<-function(variables,data=NULL,test.value,scale=FALSE,type="hist",alpha=.2){
	Value <- ..density.. <- NULL
	vars<-eval(substitute(variables),data, parent.frame())
	if(length(dim(vars))<1.5)
		vars<-d(vars)
	data <- vars
	n<-nrow(data)
	var.ind <- 1:ncol(data)
	dat<-data.frame(rep(NA,length(var.ind)*n),rep(NA,length(var.ind)*n))
	names(dat)<-c("Value","Variable")
	cnt<-1
	for(i in (1:length(var.ind))*n-(n)){
		dat[(1:n)+i,1]<-as.numeric(data[[var.ind[cnt]]])
		if(missing("test.value"))
			center <-mean(dat[(1:n)+i,1],na.rm=TRUE)
		else
			center <- test.value
		if(scale)
			dat[(1:n)+i,1]<-(dat[(1:n)+i,1]-center)/sd(dat[(1:n)+i,1],na.rm=TRUE)
		dat[(1:n)+i,2]<-names(data)[var.ind[cnt]]
		cnt<-cnt+1
	}
	dat<-na.omit(dat)
	max.diff<-max(-min(dat$Value),max(dat$Value))
	p<-ggplot(dat,aes(x=Value))
	
	hist<-geom_histogram(aes(x=Value,y=..density..))
	
	jit<-geom_jitter(aes(y=Value,x=0),alpha = alpha, colour = "#326432",  
			position = position_jitter(height = 0)) 
	box<-geom_boxplot(aes(y=Value,x=0),fill = NA, outlier.colour = "white", 
			outlier.size = 0, colour = "black", size = 0.75)
	o <- theme(axis.text.y = element_blank(), axis.title.y = element_blank())
	
	if(scale){
		p<-p+facet_grid(Variable~.)
		if(type=="hist" && !missing("test.value"))
			p <- p + hist + geom_vline(xintercept = 0,colour="red",size=1.5) +
					xlim(-max.diff,max.diff)
		else if(type=="box")
			p<-p+ jit+box+ylim(-max.diff,max.diff)+coord_flip()
		p<-p+ylab("Scaled Values")+o
	}else{
		if(ncol(vars)>1.5)
			p<-p+facet_wrap(~Variable,scales="free")+o
		else
			p<-p+ if(type=="hist") xlab(names(vars)[1]) else xlab("Value")

		if(type=="hist"){
			p<-p+hist
			if(!missing("test.value"))
				p<-p+ geom_vline(xintercept = test.value,colour="red",size=1.5)
		}else if(type=="box" && !missing("test.value"))
			p<-p+ geom_hline(yintercept = test.value,colour="red",size=1.5)+
					jit+box+o
		else if(type=="box")
			p<-p+jit+box+o + if(ncol(vars)<1.5) ylab(names(vars)[1]) else ylab("Value")
			
	}
	p
}
