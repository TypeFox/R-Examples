cvtable.plsRglm <- function(x,verbose=TRUE,...)
{
  MClassed=FALSE
  if("CV_MissClassed" %in% colnames(x[[1]])){
    if(verbose){cat("\nCV MissClassed criterion:")}
    MClassed=TRUE
    mincvMC<-function(lll){return(which.min(lll[-1,4]))} 
    mincvMCobs<-sapply(x,mincvMC)
    rescvMC<-table(factor(mincvMCobs,levels=1:max(mincvMCobs))) 
    if(verbose){print(rescvMC)}
  }
  
  if("Q2Chisq_Y" %in% colnames(x[[1]])){
    if(verbose){cat("\nCV Q2Chi2 criterion:")}
    mincvQ2Chisq<-function(lll){ if(all(lll[-1,5+2*MClassed]>lll[-1,4+2*MClassed])){return(length(lll[-1,5+2*MClassed]))} else { return(which.max(lll[-1,5+2*MClassed]<lll[-1,4+2*MClassed])-1)}}  
    mincvQ2Chisqobs<-sapply(x,mincvQ2Chisq)
    rescvQ2Chisq<-table(factor(mincvQ2Chisqobs,levels=0:max(mincvQ2Chisqobs)))
    if(verbose){print(rescvQ2Chisq)}    
    }
  
  if("PREChi2_Pearson_Y" %in% colnames(x[[1]])){
    if(verbose){cat("\nCV PreChi2 criterion:")}
    mincvPreChi2<-function(lll){return(which.min(lll[-1,6+2*MClassed]))} 
    mincvPreChi2obs<-sapply(x,mincvPreChi2)     
    rescvPreChi2<-table(factor(mincvPreChi2obs,levels=1:max(mincvPreChi2obs)))
    if(verbose){print(rescvPreChi2)}
  }
  
  if(MClassed){
    res=list(CVMC=rescvMC,CVQ2Chi2=rescvQ2Chisq,CVPreChi2=rescvPreChi2)
    class(res) <- "table.summary.cv.plsRglmmodel"
    invisible(res)    
  } else {
    res=list(CVQ2Chi2=rescvQ2Chisq,CVPreChi2=rescvPreChi2)  
    class(res) <- "table.summary.cv.plsRglmmodel"
    invisible(res)    
  }  
}

