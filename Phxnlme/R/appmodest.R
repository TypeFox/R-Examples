
"appmodest"<-function(path=NULL,model.file,out.file="out0001.txt",...)
{
  
  ctlfl<-paste(path,model.file,sep="/")
  outfl<-paste(path,out.file,sep="/")
  
  if(file.exists(ctlfl) & file.exists(outfl)){
    ctltxt<-readLines(ctlfl)
    dmptxt<-getmodest(outfl,dmpout=FALSE)
    
    ranefyes<-any(regexpr("ranef",ctltxt,fixed=TRUE)>0)
    whichef<-"fixed"
    parmtag<-"fixed"
    if(ranefyes){
      whichef<-c(whichef,"random")
      parmtag<-c(parmtag,"omega5")
    }
    whichef<-c(whichef,"error")
    parmtag<-c(parmtag,"fixed")
    
    ####4.21 SL 
    whichef<-c(whichef,"stparm")
    whichef<-c(whichef,"observe")
    
    finalest<-list(fixed=NULL,random=NULL,error=NULL)
    ########################
    
    for(i in 1:length(whichef)){
      
      ####4.21 SL 
      if(any(c(whichef[i]=="fixed",whichef[i]=="random",whichef[i]=="error"))){
        
        estmat<-dmptxt[parmtag[i]][[1]]
        if(is.vector(estmat)) estmat<-t(as.matrix(estmat))
        mathead<-colnames(estmat)
        mathead<-data.frame(start=1:length(mathead),parm=mathead,stringsAsFactors=FALSE)
      }
      
      if(whichef[i]=="fixed"){
        strtag<-"fixef("
        if(ranefyes) endtag<-"ranef(" else endtag<-"error("
      }
      if(whichef[i]=="random"){
        strtag<-"ranef("
        endtag<-"stparm("                ####4.21 SL 
      }
      if(whichef[i]=="error"){
        strtag<-"error("
        endtag<-"observe("              ####4.21 SL 
      }     
      
      ####4.21 SL 
      if(whichef[i]=="stparm"){
        strtag<-"stparm("
        endtag<-"error("
      }
      if(whichef[i]=="observe"){
        strtag<-"observe("
        endtag<-"}"
      }
      ###################
      
      modef<-gsub(" ","",gsub("\t","",ctltxt,fixed=TRUE))
      
      strln<-which(regexpr(strtag,modef,fixed=TRUE)>0)
      endln<-which(regexpr(endtag,modef,fixed=TRUE)>0)-1
      
      modef<-gsub(" ","",modef[strln:endln])
      
      ####4.21 SL 
      if(whichef[i]=="stparm"|whichef[i]=="observe")
      {if(whichef[i]=="stparm") stparm.eqt<-modef
       if(whichef[i]=="observe") observe.eqt<-modef
      }else{
        modef<-sub(strtag,"",modef,fixed=TRUE)
        indx<-length(modef)
        modef[indx]<-sub(")$","",modef[indx])
        indx<-which(regexpr("#",modef,fixed=TRUE)>0)
        modef[indx]<-substr(modef[indx],1,regexpr("#",modef[indx],fixed=TRUE)-1)
        ###SL, 4/22/15
        
        if(whichef[i]=="random") modef<-modef[nchar(modef)>1] else
          modef<-modef[nchar(modef)>2]
        
        modef<-sub("enable=","enable~",modef,fixed=TRUE)
        
        indx<-which(regexpr("=",modef,fixed=TRUE)>0)
        if(length(indx)>0){
          init<-list2mat(strsplit(modef,split="="))
          if(ncol(init)>1) init[is.na(init[,2]),2]<-""
          modef[indx]<-substr(modef[indx],1,regexpr("=",modef[indx])-1)
        }else init<-cbind(modef,"")
        indx<-which(regexpr("diag(",modef,fixed=TRUE)<0 &
                      regexpr("block(",modef,fixed=TRUE)<0 &
                      regexpr("same(",modef,fixed=TRUE)<0 )
        modef[indx]<-paste("(",modef[indx],sep="")
        modef<-gsub("[)]","",modef)
        modef<-list2mat(strsplit(modef,split="[(]"))
        if(ncol(modef)==3) modef[is.na(modef[,3]),3]<-"" else modef<-cbind(modef,"")
        modef<-matrix(modef[,c(1,3,2)],ncol=3)
        
        ##############4/23 SL
        if(whichef[i]=="random"){
          parm=modef[modef[,1]!="",3]
          parm.list=paste(parm, collapse=",")
          modef<-matrix(modef[,1:2],ncol=2)
          parm<-cbind(modef, t(list2mat(strsplit(parm.list, split=","))))
          init[init[,2]=="",2]<-init[init[,2]=="",1]
          init<-init[,2]
          init<-gsub("c", "",init)
          init<-gsub("[(]", "",init)
          init<-gsub("[)]", "",init)
          parm<-cbind(parm,init=matrix(init))
        }else{
          ################
          parm<-list2mat(strsplit(modef[,3],split=","))
          modef<-matrix(modef[,1:2],ncol=2)
          parm<-cbind(modef,parm)
          parm<-matrix(parm[,1:3],ncol=3)
          init<-gsub("c", "",init)          ###############4/28 SL
          init<-gsub("[(]", "",init)        ###############4/28 SL
          init<-gsub("[)]", " ",init)       ####<-- "[space]"
          parm<-cbind(parm,init=init[,2])}
        ###############################
        
        colnames(parm)<-c("str","fix","parm","init")
        parm<-as.data.frame(parm,stringsAsFactors=FALSE)
        parm$fix[parm$fix!=""]<-paste("(",parm$fix[parm$fix!=""],")",sep="")
        
        parm<-merge(parm,mathead,by="parm",all.x=TRUE)
        parm<-parm[order(parm$start),]
        
        #parm$stop<-vshift(matrix(parm$start,ncol=1),-1)
        #parm$stop<-parm$stop-1
        #parm$stop[nrow(parm)]<-nrow(mathead)
        
        if(whichef[i]=="random"){
          for(j in 1:nrow(parm)) parm$block.start[j]<-length(strsplit(parm$init[j],",")[[1]])
          
          parm$block.start<-parm$start-parm$block.start+1
          #parm$block.start[nrow(parm)]<-nrow(mathead)
          parm$parmlst=""
          parm$parmval=""
          for(j in 1:nrow(parm)){
            parm$parmlst[j]<-mathead[parm$start[j],2]
            if(parm$str[j]=="same")  parm$parmval[j]<-""else
              parm$parmval[j]<-paste(estmat[parm$start[j],parm$block.start[j]:parm$start[j]],collapse=",") 
          }
          #parm$parmval[regexpr(",",parm$parmval,fixed=TRUE)>0]<-paste("c(",parm$parmval[regexpr(",",parm$parmval,fixed=TRUE)>0],")",sep="")
          finalest$random<-parm
        }else{
          estmat<-t(estmat)
          estmat<-cbind(rownames(estmat),estmat)
          if(is.vector(estmat)) estmat<-t(matrix(estmat))
          colnames(estmat)<-c("parm","parmval")
          parm$parmlst<-parm$parm
          parm<-merge(parm,estmat,by="parm",all.x=TRUE)
          parm<-parm[order(parm$start),]
          parm$fix<-sub("enable~","enable=",parm$fix,fixed=TRUE)
          if(whichef[i]=="fixed") finalest$fixed<-parm 
          if(whichef[i]=="error") finalest$error<-parm
        }
        
      }  
    }
    #  ctlfl<-paste(model,"test.mdl",sep="/")
    ctltxt<-readLines(ctlfl)
    cutln<-which(regexpr("fixef",ctltxt,fixed=TRUE)>0)
    ctltxt<-ctltxt[1:(cutln-1)]
    
    for(i in 1:length(finalest)){
      if(!is.null(finalest[[i]])){
        strln<-c("\tfixef(","\tranef(","\terror(")
        endln<-c("\t)","\t)","}\n")
        
        ctlrec<-finalest[[i]]
        ctlrec$ctltxt<-""
        
        fixest<-which(ctlrec$fix %in% c("(freeze)","(enable=1)"))
        adjest<-which(!(ctlrec$fix %in% c("(freeze)","(enable=1)")))
        samestr<-which(ctlrec$str=="same")
        
        if(names(finalest)[i]=="random"){
          ##SL 4/23/15
          elm.start<-which(ctlrec$str!="")
          elm.end<-elm.start[-1]-1
          elm.end<-c(elm.end, nrow(ctlrec))
          for(j in 1:length(elm.start)){
            elm.list=paste(ctlrec$parmlst[elm.start[j]:elm.end[j]], collapse=",")
            ctlrec$ctltxt[elm.start[j]]<-paste(ctlrec$str[elm.start[j]],"(",elm.list,")",ctlrec$fix[elm.start[j]],"=c(",sep="")
          }#########################################################
          ctlrec$ctltxt[ctlrec$str==""]<-""
          
          
          fixest<-which(ctlrec$fix %in% c("(freeze)","(enable=1)"))
          adjest<-which(!(ctlrec$fix %in% c("(freeze)","(enable=1)")))
          
          vecest<-which(regexpr("c(",ctlrec$init,fixed=TRUE)>0)
          valest<-which(regexpr("c(",ctlrec$init,fixed=TRUE)<0)
          
          ctlrec$ctltxt[fixest]<-paste(ctlrec$ctltxt[fixest],ctlrec$init[fixest],sep="")
          
          ctlrec$ctltxt[adjest]<-paste(ctlrec$ctltxt[adjest],ctlrec$parmval[adjest],sep="")
          ctlrec$ctltxt[samestr]<-paste("same(",ctlrec$parmlst[samestr],")",sep="")
          ##SL 4/23/15
          ctlrec$ctltxt[elm.end]<-paste(ctlrec$ctltxt[elm.end],")",sep="")
          ctlrec$ctltxt[-elm.end]<-paste(ctlrec$ctltxt[-elm.end],",",sep="")          
          #######################
        }else{
          
          ctlrec$ctltxt[ctlrec$str!=""]<-paste(ctlrec$str[ctlrec$str!=""],"(",ctlrec$parmlst[ctlrec$str!=""],")",sep="")
          ctlrec$ctltxt[ctlrec$str==""]<-ctlrec$parmlst[ctlrec$str==""]
          
          ctlrec$ctltxt<-paste(ctlrec$ctltxt,ctlrec$fix,"=",sep="")
          
          fixest<-which(ctlrec$fix %in% c("(freeze)","(enable=1)"))
          adjest<-which(!(ctlrec$fix %in% c("(freeze)","(enable=1)")))
          
          vecest<-which(regexpr("c(",ctlrec$init,fixed=TRUE)>0)
          valest<-which(regexpr("c(",ctlrec$init,fixed=TRUE)<0)
          
          ctlrec$ctltxt[fixest]<-paste(ctlrec$ctltxt[fixest],ctlrec$init[fixest],sep="")
          
          #estlim<-matrix(,nrow(ctlrec),3)
          for(j in adjest) {
            estlim<-unlist(strsplit(ctlrec$init[j],split=","))
            if(length(estlim)==1){
              estlim<-c("",estlim,"")
            }else{
              estlim[1]<-paste( "c(",estlim[1],",",sep="")    ###############4/28 SL
              estlim[3]<-paste(",",estlim[3],")" ,sep="")     ###############4/28 SL
            }
            ctlrec$ctltxt[j]<-paste(ctlrec$ctltxt[j],estlim[1],ctlrec$parmval[j],estlim[3],sep="")
          }
        }
      }
      finalest[[i]]<-ctlrec
      
      if(names(finalest)[i]=="fixed"){
        ctltxt<-c(ctltxt,"\tfixef(")
        ctltxt<-c(ctltxt,paste("\t\t",ctlrec$ctltxt,sep=""))
        ctltxt<-c(ctltxt,"\t)")
      }
      if(names(finalest)[i]=="random"){
        ctltxt<-c(ctltxt,"\tranef(")
        ctltxt<-c(ctltxt,paste("\t\t",ctlrec$ctltxt,sep=""))
        ctltxt<-c(ctltxt,"\t)")
      }
      if(names(finalest)[i]=="error"){
        ####4.21 SL
        ctltxt<-c(ctltxt, paste("\t",stparm.eqt, sep="")) ####4.21 SL
        ctltxt<-c(ctltxt,paste("\terror(",ctlrec$ctltxt,")",sep=""))
        ctltxt<-c(ctltxt, paste("\t", observe.eqt, sep="")) ####4.21 SL
        ctltxt<-c(ctltxt,"}\n")
        
      }
    }
    
    
    return(list(strest=finalest,ctltxt=ctltxt)) 
  }else{
    write("Error: Required file(s) missing.Check command path, model.file, and out.file.","")
    return(NULL)
  }
}
