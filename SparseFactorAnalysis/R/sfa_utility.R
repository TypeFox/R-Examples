summary.sfa<-function(object,...) {
summary.inner.sfa <- function(object,interval=.9,print.dims=NULL,topbottom=5,... ){

	x<-object
	int.lo<-(1-interval)/2
	int.hi<-1-int.lo
	dim.table<-cbind(
		rowMeans(x$dim.sparse!=0),
		apply(x$dim.sparse,1,median),
		t(apply(x$dim.sparse,1,quantile,pr=c(int.lo,int.hi)))
	)
	rownames(dim.table)<-paste("Dimension",1:nrow(dim.table))
	colnames(dim.table)<-c("Prob Non-zero", "Posterior Median", paste(round(100*int.lo,3),"%",sep=""),paste(round(100*int.hi,3),"%",sep=""))

	##Table of top/bottom two dimensions
	if(length(rownames(x$coldim1) )==0) rownames(x$coldims.all)<-rownames(x$coldim1) <-paste("column",1:nrow(x$coldim1),sep="_")
	if(length(rownames(x$coldim2) )==0) rownames(x$coldim2) <-paste("column",1:nrow(x$coldim2),sep="_")

	if(length(rownames(x$rowdim1) )==0) rownames(x$rowdims.all)<-rownames(x$rowdim1) <-paste("row",1:nrow(x$rowdim1),sep="_")
	if(length(rownames(x$rowdim2) )==0) rownames(x$rowdim2) <-paste("row",1:nrow(x$rowdim2),sep="_")
	topbottom.func<-function(x2,topbottom,int.lo,int.hi){
		tempresults<-t(apply(x2,1,quantile,c(.5,int.lo,int.hi)))
		inds<-c(
		sort(tempresults[,1],decreasing=TRUE,index.return=TRUE)$ix[1:topbottom],
		rev(sort(tempresults[,1],decreasing=FALSE,index.return=TRUE)$ix[1:topbottom])
		)
		colnames(tempresults)<-c("Median", paste(round(100*int.lo,3),"%",sep=""),paste(round(100*int.hi,3),"%",sep=""))
		tempresults[inds,]
	}
	
	col1results<-topbottom.func(x$coldim1,topbottom,int.lo,int.hi)
	col2results<-topbottom.func(x$coldim2,topbottom,int.lo,int.hi)
	row1results<-topbottom.func(x$rowdim1,topbottom,int.lo,int.hi)
	row2results<-topbottom.func(x$rowdim2,topbottom,int.lo,int.hi)

tabdim1<-cbind(col1results,rownames(col2results),col2results)
cat('Estimated Dimensions\n')
cat('------------------------------------------\n')
print(round(dim.table,3)[dim.table[,1]>0,])

cat('\n First Dimension, Columns\n')
cat('------------------------------------------\n')
print(round(topbottom.func(x$coldim1,topbottom,int.lo,int.hi),3))
cat('\n First Dimension, Rows\n')
cat('------------------------------------------\n')
print(round(topbottom.func(x$rowdim1,topbottom,int.lo,int.hi),3))

cat('\n Second Dimension, Columns\n')
cat('------------------------------------------\n')
print(round(topbottom.func(x$coldim2,topbottom,int.lo,int.hi),3))

cat('\n Second Dimension, Rows\n')
cat('------------------------------------------\n')
print(round(topbottom.func(x$rowdim2,topbottom,int.lo,int.hi),3))

extra.cols<-extra.rows<-NULL


if(length(print.dims)>0){
	cat('\n Selected Columns\n')
cat('------------------------------------------\n')
	cols.out<-data.frame(matrix(NA,nrow=2*topbottom,ncol=2*length(print.dims)))
	i.count<-1
	for(i.col in 1:length(print.dims)){
				inds<-c(
		sort(x$coldims.all[,print.dims[i.col]],decreasing=TRUE,index.return=TRUE)$ix[1:topbottom],
		rev(sort(x$coldims.all[,print.dims[i.col]],decreasing=FALSE,index.return=TRUE)$ix[1:topbottom])
		)
		cols.out[,i.count]<-rownames(x$coldims.all)[inds]
		names(cols.out)[i.count]<-paste("Dim",print.dims[i.col])
		i.count<-i.count+1
		cols.out[,i.count]<-round(x$coldims.all[inds,print.dims[i.col]],3)
		names(cols.out)[i.count]<-paste("Value")

		i.count<-i.count+1

	}
	extra.cols<-cols.out
	print(cols.out)

	cat('\n Selected Rows\n')
cat('------------------------------------------\n')
	cols.out<-data.frame(matrix(NA,nrow=2*topbottom,ncol=2*length(print.dims)))
	i.count<-1
	for(i.col in 1:length(print.dims)){
				inds<-c(
		sort(x$rowdims.all[,print.dims[i.col]],decreasing=TRUE,index.return=TRUE)$ix[1:topbottom],
		rev(sort(x$rowdims.all[,print.dims[i.col]],decreasing=FALSE,index.return=TRUE)$ix[1:topbottom])
		)
		cols.out[,i.count]<-rownames(x$rowdims.all)[inds]
		names(cols.out)[i.count]<-paste("Dim",print.dims[i.col])
		i.count<-i.count+1
		cols.out[,i.count]<-round(x$rowdims.all[inds,print.dims[i.col]],3)
		names(cols.out)[i.count]<-paste("Value")

		i.count<-i.count+1

	}
	
	print(cols.out)
	extra.rows<-cols.out
	
}
	
}
summary.inner.sfa(object,...)
}


plot.sfa<-function(x,...){
plot.sfa.inner<-function(x, type="dim",main=NULL, ylabel=NULL,xlabel=NULL,dims.scatter=c(1,2),scatter.by="row",topbottom=3,...) {

makeitload<-dlgrob
makeitload<-proto

if(type%in%c("dim","scatter")==FALSE) stop("Type must be 'dim' or 'scatter' ")
#type == dim, scatter
if(type=="dim"){
if(is.null(main)) main<-"Posterior Estimate of Dimensionality by Dimension"
if(is.null(xlabel)) xlabel="Weight"
if(is.null(ylabel)) ylabel<-"Dimension"


	num.plot<-name.plot<-m.plot<-x$dim.sparse
	colnames.fill<-paste("Dimension",1:nrow(num.plot))
	for(i in 1:nrow(num.plot)) {
		num.plot[i,]<-i
		name.plot[i,]<-colnames.fill[i]
	}

	keeps<-rowMeans(x$dim.sparse!=0)>0
	num.plot<-num.plot[keeps,]
	name.plot<-name.plot[keeps,]
	m.plot<-m.plot[keeps,]
	data.plot<-data.frame(
	as.vector(name.plot),as.vector(num.plot),as.vector(m.plot)
	)
	names(data.plot)<-c("names","number","value")
	..density..<-value<-data<-number<-names<-NULL
	


g1<-ggplot(data.plot,aes(x=number,y=value,group=number)) +geom_violin(fill = "grey50", colour = "grey50",width=.5)+coord_flip()+ geom_boxplot(width=0.2,outlier.size=0)+xlab(ylabel)+ylab(xlabel)+ggtitle(main)+ scale_x_reverse(breaks=1:max(data.plot$number))

plot(g1)

}#Close out dim plot
if(type=="scatter"){
	
if(is.null(main)) main<-paste("Dimension",dims.scatter[1],"versus Dimension",dims.scatter[2])
if(is.null(xlabel)) xlabel<-paste("Dimension",dims.scatter[1])
if(is.null(ylabel)) ylabel<-paste("Dimension",dims.scatter[2])


if(scatter.by=="row"){
x1<-x$rowdims.all[,dims.scatter[1]]
y1<-x$rowdims.all[,dims.scatter[2]]
}

if(scatter.by=="col"){
x1<-x$coldims.all[,dims.scatter[1]]
y1<-x$coldims.all[,dims.scatter[2]]
}

if(length(names(x1))==0) names(x1)<-1:length(x1)
z<-names(x1)
keeps<-unique(c(
sort(x1,decreasing=TRUE,index.return=TRUE)$ix[1:topbottom],
sort(x1,decreasing=FALSE,index.return=TRUE)$ix[1:topbottom],
sort(y1,decreasing=TRUE,index.return=TRUE)$ix[1:topbottom],
sort(y1,decreasing=FALSE,index.return=TRUE)$ix[1:topbottom]
)

)
z[-keeps]<-""

q<-qplot(x1,y1)+geom_point(aes(colour=z))+xlab(ylabel)+ylab(xlabel)+ggtitle(main)
direct.label(q)
}#Close out scatter


}



plot.sfa.inner(x,...)
}
#plots the densities of the posterior

