"getmodest"<-function(outfl="out0001.txt",dmpout=TRUE,...){
  #outfl<-paste(model,"out0001.txt",sep="/")
  
  if(file.exists(outfl)){
    outtxt<-readLines(outfl)
    
    outres<-outtxt[(which(regexpr("residuals",outtxt,fixed=TRUE)>0)+1):length(outtxt)]
    outres<-gsub("\t"," ",outres)
    outres<-gsub("[ ]{2,}"," ",outres)
    outres<-gsub("^ ","",outres)
    outres<-gsub(" $","",outres)
    outres<-gsub(" ",",",outres)
    outres<-gsub("[\"]{2}",".",outres)
    outres<-gsub("\"","",outres)
    hdres<-strsplit(outres[1],",")
    outres<-list2mat(strsplit(outres[-1],","))
    dimnames(outres)[2]<-hdres
    hdres<-unlist(hdres)
    outres<-data.frame(outres,stringsAsFactors=FALSE)
    charcol<-which(regexpr("ID",hdres,fixed=TRUE)>0 | hdres=="ObsName")
    outres[,-charcol]<-matrix(as.numeric(unlist(outres[,-charcol])),ncol=length(hdres[-charcol]),byrow=FALSE)
    
    lastln<-which(regexpr("id1",outtxt,fixed=TRUE)>0)-1
    outtxt<-outtxt[1:lastln]
    outtxt<-outtxt[substr(outtxt,1,1)!="#"]
    
    stln<-which(regexpr("[A-z]{3}",substr(outtxt,1,3))>0)
    
    parm<-substr(outtxt[stln],1,18)
    cutpos<-regexpr(" ",parm,fixed=TRUE)
    indx<-which(cutpos>0)
    parm[indx]<-substr(parm[indx],1,cutpos[indx]-1)
    
    enln<-which(outtxt=="")-1
    enln<-c(1:(length(stln)-length(enln)),enln)
    datln<-data.frame(parm=parm,stln=stln,enln=enln,stringsAsFactors=FALSE)
    
    for(i in 1:length(parm)){
      assign(parm[i],outtxt[datln$stln[i]:datln$enln[i]])
    }
    dmp<-mget(parm)
    
    ss<-which(datln$stln==datln$enln)
    dmp[ss]<-gsub(" ","",dmp[ss])
    dmp[ss]<-as.numeric(substr(dmp[ss],regexpr("=",dmp[ss],fixed=TRUE)+1,nchar(dmp[ss])))
    
    ss<-which(parm %in% c("fixed","varFixefInf"))
    for(i in ss){
      xxx<-dmp[[i]][-1]
      xxx<-sub("\t"," ",xxx)
      nmxxx<-gsub(" ","",substr(xxx,regexpr("#",xxx,fixed=TRUE)+1,nchar(xxx)))
      xxx<-as.numeric(gsub(" ","",substr(xxx,regexpr("=",xxx,fixed=TRUE)+1,regexpr("#",xxx,fixed=TRUE)-1)))
      names(xxx)<-nmxxx
      dmp[[i]]<-xxx
    }
    
    ss<-which(parm=="eigenvalues")
    if(length(ss)>0){
      dmp[[ss]]<-dmp[[ss]][-1]
      dmp[[ss]]<-gsub("\t","",dmp[[ss]])
      dmp[[ss]]<-as.numeric(gsub(" ","",dmp[[ss]]))
    }
    
    ss<-which(parm %in% c("varFix","omega5","omega5chol","stdOmega5"))
    for(i in ss){
      xxx<-dmp[[i]]
      xxx<-gsub("\t$","",xxx)
      xxx<-gsub("\t",",",xxx)
      hdxxx<-unlist(strsplit(gsub(" ","",substr(xxx[1],regexpr("#",xxx[1],fixed=TRUE)+1,nchar(xxx[1]))),","))
      xxx<-list2mat(strsplit(gsub(" ","",xxx[-1]),","))
      nhd<-length(hdxxx)
      for(ii in 1:nhd){
        for(jj in 1:nhd){
          xxx[ii,jj]<-xxx[jj,ii]
        }
      }
      xxx<-matrix(as.numeric(as.vector(xxx)),ncol=nhd,byrow=TRUE)
      colnames(xxx)<-hdxxx
      rownames(xxx)<-hdxxx
      dmp[[i]]<-xxx
    }
    
    dmp$residuals<-outres
    
    if(dmpout){
      dput(dmp,file="outdmp.txt")
    }
    return(dmp)
    
  }else{
    write(paste("Error: Output file",outfl,"does not exist."),"")
    return(NULL)
  }
}
