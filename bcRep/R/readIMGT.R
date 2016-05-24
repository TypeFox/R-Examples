## Julia Bischof
## 02-22-2016

readIMGT<-function(data,filterNoResults=TRUE){
  if(substr(strsplit(data,split="/")[[1]][length(strsplit(data,split="/")[[1]])],1,2) %in% c("3_","6_")){
    newcolnames<-gsub(" |-","_",strsplit(readLines(data)[1],split="\t")[[1]])
    newcolnames<-gsub("'","",newcolnames)
    
    if(length(strsplit(readLines(data)[2],split="\t")[[1]])==3){
      temp<-lapply(readLines(data),function(x){strsplit(x,split="\t")[[1]]})
      for(i in 1:length(temp)){
        if(length(temp[[i]])==3){
          temp[[i]]<-c(temp[[i]],rep(" ",(length(newcolnames)-3)))
        }
      }
      tab<-t(data.frame(temp[2:length(temp)], stringsAsFactors = F))
      rownames(tab)<-as.numeric(seq(1, (length(temp)-1), 1))
    }else{
      tab<-data.frame(matrix(NA, nrow=0,ncol=length(newcolnames)))
      tab<-rbind(tab,read.table(data,fill=TRUE,header=F,skip=1,sep="\t",as.is=TRUE))
    }
    colnames(tab)<-newcolnames
  }else{
    tab<-read.table(data,fill=TRUE,header=TRUE,sep="\t",as.is=TRUE,check.names=F)
    colnames(tab)<-gsub(" |-","_",gsub("'","",colnames(tab)))
  }
  tab<-tab[,-ncol(tab)]
  if(filterNoResults==T){
    out<-grep("No results",tab[,3])
    if(length(out)>0){
      tab<-tab[-out,]
      cat("...", length(out), "sequences were excluded\n")
    }else{
      cat("... no sequences were excluded\n")
    }      
  }
  return(tab)
}


