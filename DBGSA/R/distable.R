distable<-function(fd,num,setdis,Meth,resultname)
{
a<-list()
lines=readLines(fd,n=-1)
for(i in 1:length(lines)){
ar=unlist(strsplit(lines[i],"\t"))

if(is.null(a[[ar[1]]])){
	a[[ar[1]]]=list();
	a[[ar[1]]][[1]]<-as.numeric(ar[3:length(ar)])
}else{
	len=length(a[[ar[1]]])
	len=len+1
	a[[ar[1]]][[len]]<-as.numeric(ar[3:length(ar)])
}
}

b<-list()
name<-names(a)

for(i in 1:length(name)){
	rownum<-length(a[[name[i]]])
	colnum<-length(a[[name[i]]][[1]])
	mat<-matrix(nrow=rownum,ncol=colnum)
	for(j in 1:rownum){
		mat[j,]<-a[[name[i]]][[j]]
	}
	b[[name[i]]]<-mat
}

result<-matrix(nrow=length(name),ncol=2)
for(i in 1:length(name)){
	C<-b[[name[i]]]
	
result[i,1]<-length(a[[name[i]]])
result[i,2]<-setdis(C,num,Meth)

}
df<-data.frame(GOid=name,genenum=result[,1],dis0=result[,2])
write.table(df,file=resultname,row.names=FALSE,quote=FALSE)
}