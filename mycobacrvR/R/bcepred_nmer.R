bcepred_nmer <-
function(clas)
{
i1=1
i2=length(dir(pattern=".Rdata"))
 for(i in i1:i2)
  {
     load(dir(pattern=".Rdata")[i])
     temp=ls(pattern=clas)
     temp2=lapply(temp,function(x)
     {tempm=get(x);
      lapply(tempm,function(x)
      {x@sequence}
             )
     }
                  )  
     
     org=NULL;for(i in 1:length(temp)){org=append(org,strsplit(temp[i],split="_",fixed=T)[[1]][1])}
     rv=which(org=="mtbh37rv")
    temp4=NULL; for(i in 1:length(temp2[[rv]])){for(j in 1:length(temp2)){temp4=c(temp4, length(intersect(temp2[[j]],temp2[[rv]][i])))}}
     
     temp5=matrix(temp4, nrow=length(temp2[[rv]]), ncol=length(temp), byrow=T)
     colnames(temp5)=substr(temp,1,8)
     epitope_seq= as.matrix(temp2[[rv]],ncol=1,byrow=T)
     gi_number=as.matrix(strsplit(temp[rv],split="epitopes",fixed=T)[[1]][2],ncol=1,byrow=T)
     orthologs=as.matrix(length(temp)-1,ncol=1,byrow=T)
     epitope_length=NULL;for(i in 1:length(get(temp[rv])))
     {epitope_length=append(epitope_length,get(temp[rv])[[i]]@length)}
     class=as.matrix("bcepred",ncol=1,byrow=T);
     conservation_ratio=NULL;for(k in 1:nrow(temp5))
     {conservation_ratio=append(conservation_ratio,sum(temp5[k,1:ncol(temp5)])-1)}
     temp6=as.data.frame(temp5)
     temp6=cbind(gi_number,orthologs,epitope_length,epitope_seq,conservation_ratio,class,temp6)
     colnames(temp6)[1:6]=c("ginumber","Number_of_Orthologs","Epitope_Length","Epitope_Sequence","Epitope_Conservation_Ratio","Class")
     temp7=temp6[1:nrow(temp6),1:6]
     out <- data.frame(lapply(temp7, function(x) factor(unlist(x))))
     write.table(out,file=paste("bcepred_",strsplit(temp[rv],split="epitopes",fixed=T)[[1]][2],".txt",sep=""),sep="\t",row.names=F)
     write.table(temp5,file=paste("bcepred_",strsplit(temp[rv],split="epitopes",fixed=T)[[1]][2],"_metadata.txt",sep=""),sep="\t",row.names=F)
     #save.image(file=paste("bcepred_",strsplit(temp[rv],split="epitopes",fixed=T)[[1]][2],".RData",sep=""))
rm(list= ls()[!(ls() %in% c('clas'))])     
}
}
