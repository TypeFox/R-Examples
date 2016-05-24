plot.cor.matrix<-function(x,y=NULL,size=10,...){
	x.o <- y.o <- NULL
	x<-as.matrix(x)
	n.r<-nrow(x)
	n.c<-ncol(x)
	dat<-as.data.frame(matrix(NA,n.r*n.c,3))
	cnt<-1
	for(i in 1:n.r)
		for(j in 1:n.c){
			dat[cnt,1]<-x[i,j]
			dat[cnt,2]<-rownames(x)[i]
			dat[cnt,3]<-colnames(x)[j]
			cnt<-cnt+1
		}
	dat<-cbind(dat,0)
	dat<-cbind(dat,0)
	names(dat)<-c("cor","row","col","x.o","y.o")	
	p<-qplot(x.o,y.o,data=dat,facets=row~col,size=abs(cor),colour=cor)+
			scale_size(limits=c(0,1),range=c(0,size))+xlim(c(-1,1))+ylim(c(-1,1))
	p<-p+scale_colour_gradient(limits=c(-1,1),low="#B71B1A",high="#3B4FB8")
	o <- theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
			axis.ticks = element_blank(), axis.text.y = element_blank(), 
			axis.text.x = element_blank(), axis.title.y = element_blank(), 
			axis.title.x = element_blank())
	p+o
}
