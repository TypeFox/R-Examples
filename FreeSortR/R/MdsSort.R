
################################################################
#________________      Function  MdsSort        ________________
#                  Mds of free sorting data
################################################################

MdsSort<-function(Part,ndim=2,nboot=0,metric=FALSE,ties="primary",itmax=5000,eps=1e-06){
  
  if (!class(Part)=="SortingPartition"){
    return("This is not an object of class SortingPartition")
  }  
  else
  {
    nstimuli<-Part@nstimuli
    nsubjects<-Part@nsubjects
    
    ListDissim<-Dissimil(Part)
    
    MatDissim<-apply(simplify2array(ListDissim),c(1,2),'sum')
    
    #library(smacof)
    
    res<-MdsDiss(MatDissim,ndim=ndim,metric=metric,ties,itmax,eps)
    Config<-res$Config
    rownames(Config)<-Part@LabStim
    colnames(Config)<-paste(rep("Dim"),1:ndim)
    Stress<-res$Stress
    Percent<-res$Percent[1:ndim]
    
    if (nboot>0){
      
      #library(vegan)
      
      ResBoot<-vector("list",nstimuli)
      for (b in 1:nboot){
        panel<-sample(1:nsubjects,replace=TRUE)
        MatDistTotPanel= matrix(0, nstimuli,nstimuli)  
        for (k in 1:nsubjects){                       
          MatDistTotPanel<-MatDistTotPanel+ListDissim[[panel[k]]]
        }
        if (metric==TRUE){
          res<-smacofSym(MatDistTotPanel,ndim=ndim,type="ratio",ties="primary",init=Config,verbose=F,itmax=itmax,eps=eps)
        } else {
          res<-smacofSym(MatDistTotPanel,ndim=ndim,type="ordinal",ties="primary",init=Config,verbose=F,itmax=itmax,eps=eps)
        }
        ConfigP<-res$conf
        
        res<-procrustes(Config,ConfigP)
        for (i in 1:nstimuli){
          ResBoot[[i]]<-rbind(ResBoot[[i]],res$Yrot[i,])
        } 
      }
      
      return(new("SortingMds",nstimuli=nstimuli,nsubjects=nsubjects,LabStim=Part@LabStim,LabSubj=Part@LabSubj,ndim=ndim,Config=Config,Percent=Percent,Stress=Stress,ResBoot=ResBoot))
      #return(list(Config=Config,Percent=Percent,Stress=Stress,ResBoot=ResBoot))
      
    }else{
      ResBoot<-NULL
      return(new("SortingMds",nstimuli=nstimuli,nsubjects=nsubjects,LabStim=Part@LabStim,LabSubj=Part@LabSubj,ndim=ndim,Config=Config,Percent=Percent,Stress=Stress))
      #return(list(Config=Config,Percent=Percent,Stress=Stress))
    }
  }
  
}  
