MDocWordTable<-function(MDocTerm,MDocVar,Fmin=rep(5,length(MDocTerm)),Dmin=rep(5,length(MDocTerm)),
idiom=rep("en",length(MDocTerm)),stop.word.user=vector(mode="list",length(MDocTerm)),
stop.word.tm=rep(FALSE,length(MDocTerm)),num.agg=NULL,Fmax=NULL){
      
  	MTables<-list()
	DocTermVar<-list()
	MTagregSep<-list()
      gpo<-vector()
 for(i in 1:length(MDocTerm)){
      num.agg=num.agg
      DocTerm<-as.data.frame(MDocTerm[[i]])
	DocVar<-as.data.frame(MDocVar[[i]])
      rownames(DocVar)<-rownames(DocTerm)
	dtmM<-as.matrix(DocTerm)
      FreqWord<-apply(dtmM,MARGIN=2,FUN=sum)
      dtmA<-dtmM
      dtmA[dtmM[,]>0]<-1
      NumDoc<-apply(dtmA,MARGIN=2,FUN=sum)
      DocTermR<-DocTerm[,which(FreqWord>=Fmin[i]&NumDoc>=Dmin[i])]
        if(!is.null(Fmax)){
         FrecMax<-apply(DocTermR,2,sum)
         PalSel<-which(FrecMax>=Fmax)
        DocTermR<-DocTermR[,-PalSel]
       }

        if (!is.null(stop.word.user[[i]])){
               WordUser<- as.vector(stop.word.user[[i]])
                DocTermR<-DocTermR[,which(!colnames(DocTermR)%in%WordUser)]
        }
        if(stop.word.tm[i]){
                stopword<-stopwords(idiom[i])
                DocTermR<-DocTermR[,which(!colnames(DocTermR)%in%stopword)]
        }
        dim(DocTermR)       
        DocTermR<-DocTermR[apply(DocTermR,1,sum)>0,]
        Dcol<-ncol(DocTermR)
        gpo[i]<-Dcol
        MTables[[i]]<-DocTermR
        DocVarR<-as.data.frame(DocVar[which(rownames(DocVar)%in%rownames(DocTermR)),])
        dim(DocVarR)
	 if(length(DocVarR)==1)
        colnames(DocVarR)<-colnames(DocVar)
        DocTermV<-cbind.data.frame(DocTermR,DocVarR)
        DocTermVar[[i]]<-DocTermV   
        base<-DocVarR
     if(!is.null(num.agg)){
        AggregVar<-num.agg
     	if(length(AggregVar)==1)
		AggregVar<-AggregVar
   	if(length(AggregVar)>1)
		AggregVar<-AggregVar[i]
    	if(is.character(AggregVar))
           AggregVar<-which(colnames(base)%in%AggregVar)
          if(is.numeric(AggregVar))
             AggregVar<-AggregVar 
         agg<-(base[,AggregVar])
         DocTermRA<-DocTermR
         dis.X<-tab.disjonctif(agg)
         Tagreg<-t(DocTermRA)%*%dis.X
         Tagreg<-t(Tagreg)
		if(i==1) Mtable=Tagreg
             else Mtable<-cbind.data.frame(Mtable,Tagreg)
          MTagregSep[[i]]<-Tagreg
     
 	}
      else{
               if(i==1) Mtable=DocTermR[,apply(DocTermR,2,sum)>0]
                
             else{   
			 DocTermR<-DocTermR[which(rownames(DocTermR)%in%rownames(Mtable)),]
			 Mtable<-Mtable[which(rownames(Mtable)%in%rownames(DocTermR)),]
                   Mtable<-Mtable[(apply(Mtable,1,sum)>0)&(apply(DocTermR,1,sum)>0),]
                   Mtable<-Mtable[,apply(Mtable,2,sum)>0]  
                   DocTermR<-DocTermR[(apply(Mtable,1,sum)>0)&(apply(DocTermR,1,sum)>0),]
                   DocTermR<-DocTermR[,apply(DocTermR,2,sum)>0] 
			 Mtable<-cbind.data.frame(Mtable,DocTermR)
                   gpo[i]<-ncol(DocTermR)
            
          }  
		MTagregSep=NULL
       }
  }  
          
    res<-list(MDocWord=Mtable, ncolTs=gpo,MTagregSep=MTagregSep)
return(res)
}

