save.var.txt <- function(wave.var.list)
{

if(class(wave.var.list)!="Wave Variance")
          stop("The type of the data is not correct, please check!")

n.levels<-length(wave.var.list)/3
n.regions<-length(wave.var.list[[1]])


# To test if wave.var.list is of version 1 or 2 and after

list.attr<-names(attributes(wave.var.list))
version<-1

for(i in 1:length(list.attr)){
if(list.attr[i]=="version") version<-attr(wave.var.list, "version")
}

method <- attr(wave.var.list, "method")
wf <- attr(wave.var.list, "wavelet")
boundary <- attr(wave.var.list, "boundary")
if(version!=1) proc.length<- attr(wave.var.list, "proc.length")

for(i in 1:n.levels){
name.txt<-paste("wave_var_mat_level_",i,sep="")
name.txt<-paste(name.txt,".txt",sep="")

if(version!=1) param.txt<-c(method,wf,boundary,n.levels,proc.length,n.regions,"data begins line 4")
if(version==1) param.txt<-c(method,wf,boundary,n.levels,"data begins line 4")

write.table(version,name.txt,row.names=FALSE,col.names=FALSE,eol=" ",quote=FALSE)
write.table(" ",name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol="\n",quote=FALSE)
write.table(class(wave.var.list),name.txt,row.names=FALSE,col.names=FALSE,eol=" ",quote=FALSE,append=TRUE)
write.table(" ",name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol="\n",quote=FALSE)
write.table(param.txt,name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol=" ",quote=FALSE)
write.table(" ",name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol="\n",quote=FALSE)
write.table(wave.var.list[[i]],name.txt,row.names=FALSE,col.names=FALSE,append=TRUE)

name.txt<-paste("wave_var_lower_mat_level_",i,sep="")
name.txt<-paste(name.txt,".txt",sep="")

write.table(version,name.txt,row.names=FALSE,col.names=FALSE,eol=" ",quote=FALSE)
write.table(" ",name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol="\n",quote=FALSE)
write.table("Wave Variance lower bound",name.txt,row.names=FALSE,col.names=FALSE,eol=" ",quote=FALSE,append=TRUE)
write.table(" ",name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol="\n",quote=FALSE)
write.table(param.txt,name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol=" ",quote=FALSE)
write.table(" ",name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol="\n",quote=FALSE)
write.table(wave.var.list[[i+4]],name.txt,row.names=FALSE,col.names=FALSE,append=TRUE)

name.txt<-paste("wave_var_upper_mat_level_",i,sep="")
name.txt<-paste(name.txt,".txt",sep="")

write.table(version,name.txt,row.names=FALSE,col.names=FALSE,eol=" ",quote=FALSE)
write.table(" ",name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol="\n",quote=FALSE)
write.table("Wave Variance upper bound",name.txt,row.names=FALSE,col.names=FALSE,eol=" ",quote=FALSE,append=TRUE)
write.table(" ",name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol="\n",quote=FALSE)
write.table(param.txt,name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol=" ",quote=FALSE)
write.table(" ",name.txt,row.names=FALSE,col.names=FALSE,append=TRUE,eol="\n",quote=FALSE)
write.table(wave.var.list[[i+8]],name.txt,row.names=FALSE,col.names=FALSE,append=TRUE)

}

}

