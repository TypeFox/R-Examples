cvtable.plsR <- function(x,verbose=TRUE,...)
{
  MClassed=FALSE
  if("CV_MissClassed" %in% colnames(x[[1]])){
    if(verbose){cat("\nCV MissClassed criterion:")}
    MClassed=TRUE
    mincvMC<-function(lll){return(which.min(lll[-1,3]))} 
    mincvMCobs<-sapply(x,mincvMC)
    rescvMC<-table(factor(mincvMCobs,levels=1:max(mincvMCobs))) 
    if(verbose){print(rescvMC)}
  }
  
  if("Q2_Y" %in% colnames(x[[1]])){
    if(verbose){cat("\nCV Q2 criterion:")}
    mincvQ2<-function(lll){ if(all(lll[-1,4+2*MClassed]>lll[-1,3+2*MClassed])){return(length(lll[-1,4+2*MClassed]))} else { return(which.max(lll[-1,4+2*MClassed]<lll[-1,3+2*MClassed])-1)}}  
    mincvQ2obs<-sapply(x,mincvQ2)
    rescvQ2<-table(factor(mincvQ2obs,levels=0:max(mincvQ2obs)))   
    if(verbose){print(rescvQ2)}    
  }
  
  if("PRESS_Y" %in% colnames(x[[1]])){
    if(verbose){cat("\nCV Press criterion:")}
    mincvPress<-function(lll){return(which.min(lll[-1,5+2*MClassed]))} 
    mincvPressobs<-sapply(x,mincvPress)     
    rescvPress<-table(factor(mincvPressobs,levels=1:max(mincvPressobs)))
    if(verbose){print(rescvPress)}
  }
  
  if(MClassed){
    res=list(CVMC=rescvMC,CVQ2=rescvQ2,CVPress=rescvPress)
    class(res) <- "table.summary.cv.plsRmodel"
    invisible(res)    
  } else {
    res=list(CVQ2=rescvQ2,CVPress=rescvPress)  
    class(res) <- "table.summary.cv.plsRmodel"
    invisible(res)    
  }  
}

