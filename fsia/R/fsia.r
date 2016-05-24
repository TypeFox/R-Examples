# conc numeric or character
# key numeric, rownumber, or character, key
read.formscanner<-function(file,conc=NULL,keyline=NULL,keyfile=NULL,keydata=NULL,id)
{
	data<-read.csv2(file)
	key<-NULL
	if (is.numeric(conc)) {
		tmp<-c()
		for (i in conc) tmp<-paste(tmp,data[,i],sep="")
		data[,conc[1]]<-tmp
		for (i in 2:length(conc)) data[,conc[i]]<-c()
	}
	if (!is.numeric(conc)) {
		tmp<-c()
		for (i in conc) tmp<-paste(tmp,data[,i],sep="")
		data[,conc[1]]<-tmp
		for (i in 2:length(conc)) data[,conc[i]]<-c()
	}
	colnames(data)[colnames(data)==id]<-"id"
	if (!is.null(keyline)) {
		key<-data[keyline,]
		data<-data[-keyline,]
	}
	if (!is.null(keyfile)) key<-read.csv2(keyfile)
	if (!is.null(keydata)) key<-keydata
	#if (!is.null(delete)) data[,colnames(data)%in%delete]<-NULL
	sel<-colnames(key)[colnames(key)%in%colnames(data)]
	if (!is.null(key)) key<-key[,sel]
	out<-list(data=data,key=key)
	class(out)<-"fsdata"
	return(out)
}


resp2binary<-function(obj,columns)
{
	key<-obj$key
	if (is.null(key)) stop("key is required")
	if (is.numeric(columns)) item<-colnames(obj$data)[columns]
	else item<-columns
	data<-obj$data
	out<-matrix(NA,nrow(data),length(columns))
	for (i in 1:length(columns)) {
		out[,i]<-(as.character(data[,item[i]])==as.character(key[,item[i]]))*1
	}
	data[,columns]<-out
	return(data)
}


freq<-function(obj,columns,perc=FALSE)
{
	if (is.numeric(columns)) item<-colnames(obj$data)[columns]
	else item<-columns
	if (!is.null(obj$key)) key<-as.matrix(obj$key[,item])
	else key<-obj$key
	out<-list()
	j<-1
	for (i in columns) {
		tab<-table(obj$data[,i])
		if (perc) {
			tab<-tab/nrow(obj$data)*100
			tab<-round(tab,2)
			#tab<-paste(tab,"%",sep="")
		}
		out[[j]]<-list(item=item[j],tab=tab,key=key[j])
		j<-j+1
	}
	class(out)<-"frlist"
	return(out)
}


print.frlist<-function(x, ...)
{
	for (i in 1:length(x)) {
		cat("\n============== ")
		cat(x[[i]]$item)
		cat(" ==============\n")
		tab<-x[[i]]$tab
		names(tab)[names(tab)==x[[i]]$key]<-paste(names(tab)[names(tab)==x[[i]]$key],"*",sep="")
		print(tab)
	}
	cat("\n")
}

plot.frlist<-function(x, display=TRUE, ...)
{
	devAskNewPage(ask = TRUE)
	for (i in 1:length(x)) {
		tab<-x[[i]]$tab
		colour<-rep(2,dim(tab))
		if (!is.null(x[[i]]$key)) colour[names(tab)==x[[i]]$key]<-3
		bp<-barplot(tab,col=colour,main=x[[i]]$item,ylim=c(0,max(tab)*1.2))
		if (display) {
			text(x=bp,y=(tab+max(tab)*0.02),labels=tab,adj = c(0.5, 0))
		}
	}
	devAskNewPage(ask = FALSE)
}




person.stat<-function(obj,columns)
{
	data01<-resp2binary(obj=obj,columns=columns)
	score<-rowSums(data01[,columns])
	out<-data.frame(id=obj$data$id,score=score,count=length(columns),perc=(score/length(columns)*100))
	return(out)
}



item.stat<-function(obj,columns)
{
	data01<-resp2binary(obj=obj,columns=columns)
	score<-colSums(data01[,columns])
	out<-data.frame(item=names(score),score=score,count=nrow(obj$data),perc=(score/nrow(obj$data)*100))
	rownames(out)<-NULL
	return(out)
}



report<-function(obj,columns,whichid,grid=TRUE,main="",las=0)
{
	if (is.numeric(columns)) item<-colnames(obj$data)[columns]
	else item<-columns
	n<-length(columns)
	h<-rep(0,(length(whichid)+2))
	names(h)<-c("item",whichid,"key")
	bp<-barplot(h,axes=FALSE,xlab="",ylab="", border = NA ,ylim=c(0,n),main=main,las=las)
	if (!is.null(obj$key)) key<-as.matrix(obj$key[,item])
	else key<-obj$key
	resp<-as.matrix(obj$data[obj$data$id%in%whichid,columns])
	text(bp[1],n:1-0.5,item)
	for (i in seq_along(whichid)) {
		colour<-rep(2,n)
		colour[resp[i,]==key]<-3
		text(bp[i+1],n:1-0.5,resp[i,],col=colour)
	}
	text(bp[i+2],n:1-0.5,key)
	if (grid) abline(h=(0:n))
}



