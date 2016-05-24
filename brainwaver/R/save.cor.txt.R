save.cor.txt <- function(wave.cor.list)
{

if(class(wave.cor.list)!="Wave Correlation")
          stop("The type of the data is not correct, please check!")

n.levels<-length(wave.cor.list)/3
n.regions<-dim(wave.cor.list[[1]])[1]

# To test if wave.cor.list is of version 1 or 2 and after

list.attr<-names(attributes(wave.cor.list))
version<-1

for(i in 1:length(list.attr)){
if(list.attr[i]=="version") version<-attr(wave.cor.list, "version")
}

method <- attr(wave.cor.list, "method")
wf <- attr(wave.cor.list, "wavelet")
boundary <- attr(wave.cor.list, "boundary")
if(version!=1) proc.length<- attr(wave.cor.list, "proc.length")


for(i in 1:n.levels){
name.txt<-paste("wave_cor_mat_level_",i,sep="")
name.txt<-paste(name.txt,".txt",sep="")

if(version!=1) param.txt<-c(method,wf,boundary,n.levels,proc.length,n.regions,"data begins line 4")
if(version==1) param.txt<-c(method,wf,boundary,n.levels,"data begins line 4")

write.table(version,name.txt,row.names=FALSE,col.names=FALSE,eol=" ",quote=FALSE)
write.table(" ",name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol="\n",quote=FALSE)
write.table(class(wave.cor.list),name.txt,row.names=FALSE,col.names=FALSE,eol=" ",quote=FALSE,append=TRUE)
write.table(" ",name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol="\n",quote=FALSE)
write.table(param.txt,name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol=" ",quote=FALSE)
write.table(" ",name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol="\n",quote=FALSE)
write.table(wave.cor.list[[i]],name.txt,row.names=FALSE,col.names=FALSE,append=TRUE)

name.txt<-paste("wave_cor_lower_mat_level_",i,sep="")
name.txt<-paste(name.txt,".txt",sep="")

write.table(version,name.txt,row.names=FALSE,col.names=FALSE,eol=" ",quote=FALSE)
write.table(" ",name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol="\n",quote=FALSE)
write.table("Wave Correlation lower bound",name.txt,row.names=FALSE,col.names=FALSE,eol=" ",quote=FALSE,append=TRUE)
write.table(" ",name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol="\n",quote=FALSE)
write.table(param.txt,name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol=" ",quote=FALSE)
write.table(" ",name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol="\n",quote=FALSE)
write.table(wave.cor.list[[i+4]],name.txt,row.names=FALSE,col.names=FALSE,append=TRUE)

name.txt<-paste("wave_cor_upper_mat_level_",i,sep="")
name.txt<-paste(name.txt,".txt",sep="")

write.table(version,name.txt,row.names=FALSE,col.names=FALSE,eol=" ",quote=FALSE)
write.table(" ",name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol="\n",quote=FALSE)
write.table("Wave Correlation upper bound",name.txt,row.names=FALSE,col.names=FALSE,eol=" ",quote=FALSE,append=TRUE)
write.table(" ",name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol="\n",quote=FALSE)
write.table(param.txt,name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol=" ",quote=FALSE)
write.table(" ",name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol="\n",quote=FALSE)
write.table(wave.cor.list[[i+8]],name.txt,row.names=FALSE,col.names=FALSE,append=TRUE)

}

}

