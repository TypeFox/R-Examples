"savematrix" <-
structure(function(n,filename,twomode=1,cont=FALSE){
if(length(grep(patt="w32",x=version["os"]))){
	eol<-"\n"
}else{eol<-"\r\n"}
if ((dim(n)[1] == dim(n)[2]) & (twomode!=2))
{ 
  verNames<-rownames(n)
  if(is.null(verNames))verNames<-1:dim(n)[1]
  verNamesTable<-table(verNames)
  if(max(verNamesTable)>1){
  	duplicateName<-names(which(verNamesTable>1))
  	for(i in duplicateName){
  		verNames[verNames==i]<-paste(i,1:verNamesTable[i],sep="")
  	}
  }
  cat(paste("*Vertices",dim(n)[1]),eol, file = filename,append=cont);
  cat(paste(seq(1,length=dim(n)[1]),' "',verNames,'"',eol,sep=""), file = filename,append=TRUE);
  cat("*Matrix",eol, file = filename,append=TRUE);
  write.table(n,file=filename,eol=eol,row.names = FALSE, col.names = FALSE,append=TRUE)
}else
{ 
  verRowNames<-rownames(n)
  if(is.null(verRowNames))verRowNames<-1:dim(n)[1]
  verRowNamesTable<-table(verRowNames)
  if(max(verRowNamesTable)>1){
  	duplicateRowName<-names(which(verRowNamesTable>1))
  	for(i in duplicateRowName){
  		verRowNames[verRowNames==i]<-paste(i,1:verRowNamesTable[i],sep="")
  	}
  }
  verColNames<-colnames(n)
  if(is.null(verColNames))verColNames<-1:dim(n)[2]
  verColNamesTable<-table(verColNames)
  if(max(verColNamesTable)>1){
  	duplicateColName<-names(which(verColNamesTable>1))
  	for(i in duplicateColName){
  		verColNames[verColNames==i]<-paste(i,1:verColNamesTable[i],sep="")
  	}
  }
  cat(paste("*Vertices",sum(dim(n)),dim(n)[1]),eol, file = filename,append=cont);
  cat(paste(1:dim(n)[1],' "',verRowNames,'"',eol,sep=""), file = filename,append=TRUE);
  cat(paste(seq(dim(n)[1]+1,length=dim(n)[2]),' "',verColNames,'"',eol,sep=""), file = filename,append=TRUE);
  cat("*Matrix",eol, file = filename, append=TRUE);
  write.table(n,file=filename,eol=eol,row.names = FALSE, col.names = FALSE,append=TRUE)
} }
, comment = "Save matrix to file that can be read by Pajek (as *Matrix)")
