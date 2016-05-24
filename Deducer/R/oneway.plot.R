oneway.plot<-function(formula,data=NULL,alpha=.2,
		box=TRUE,points=TRUE,scale=FALSE){
	Factor <- Value <- NULL
	n<-nrow(data)
	variables<-eval(formula[[2]],data,parent.frame())
	factor.var<-eval(formula[[3]],data,parent.frame())
	if(length(dim(variables))<1.5){
		variables<-d(variables)
		fn<-formula[[2]]
		names(variables)<-if(is.call(fn)) format(fn) else as.character(fn)
	}else
		variables<-as.data.frame(variables)
	if(length(dim(factor.var))>1.5)
		factor.var <- factor.var[,1]
	n.var<-ncol(variables)
	dat<-data.frame(rep(NA,n.var*n),rep(NA,n.var*n),
			rep(NA,n.var*n))
	names(dat)<-c("Value","Factor","Variable")
	cnt<-1
	for(i in (1:n.var)*n-(n)){
		dat[(1:n)+i,1]<-as.numeric(variables[[cnt]])
		if(scale)
			dat[(1:n)+i,1]<-(dat[(1:n)+i,1]-mean(dat[(1:n)+i,1],na.rm=TRUE))/
					sd(dat[(1:n)+i,1],na.rm=TRUE)
		dat[(1:n)+i,3]<-names(variables)[cnt]
		cnt<-cnt+1
	}
	dat[,2]<-rep(factor.var,n.var)
	dat[,2]<-as.factor(dat[,2])
	dat[,3]<-as.factor(dat[,3])
	dat<-na.omit(dat)
	p<-ggplot(dat, aes(Factor,Value))
	if(points){
		p<-p+geom_jitter(aes(colour=Factor),alpha = alpha,
				position=position_jitter(height=0))+theme(legend.position = "none")
		if(length(levels(dat$Factor))>2 & length(levels(dat$Factor))<9)
			p<-p+scale_colour_brewer(palette="Dark2") 
		if(box)
			p<-p+geom_boxplot(fill = NA, outlier.colour="white", outlier.size=0,
					colour="black",size=.75)
	}else
	if(box)
		p<-p+geom_boxplot(colour="black",size=.75)
	else
		stop("Either box or points must be TRUE") 
	if(!scale){
		if(n.var==1)
			p<-p+ ylab(dat[1,3]) 
		else 
			p<-p+facet_wrap(~Variable,scales="free_y")+ylab("")
	}else{
		if(n.var==1) 
			p<-p+ylab(paste("Scaled",dat[1,3])) 
		else 
			p<-p+facet_grid(Variable~.,scales="free_y")+ylab("Standardized Value")
	}
	fn<-formula[[3]]
	fact.name<-if(is.call(fn)) format(fn) else as.character(fn)
	p<-p+xlab(fact.name)#+coord_flip()
	p
}