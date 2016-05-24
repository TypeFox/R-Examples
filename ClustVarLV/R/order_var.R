order_var <-
function(ordr,listind) 
{
  ordr <- c(ordr,listind);
  dupli<-duplicated(ordr,fromLast=FALSE)
  duplibis<-duplicated(ordr,fromLast=TRUE)
  if (sum(dupli)>0) {
    idem<-which(dupli)
    dupli[(length(ordr)-length(listind)+1):(length(ordr))]<-as.logical(1-dupli[(length(ordr)-length(listind)+1):(length(ordr))])
    dif<-which(dupli)    
    if(length(dif)>0) {
      num<-ordr[dif]
      ou<-which(duplibis)     
      tempo=matrix(0,nrow=1,ncol=length(ordr)-length(idem))
      tempo[1:max(ou)]<-ordr[1:max(ou)]
      tempo[max(ou)+1]<-num
      if (length(tempo)>(max(ou)+1))  tempo[(max(ou)+2):(length(ordr)-length(idem))]<-ordr[(max(ou)+1):(length(ordr)-length(listind))]  
      ordr<-tempo
    }
    else {                             
      ordr<-ordr[-((length(ordr)-length(listind)+1):length(ordr))]                                    
    }     
  }
  
return(ordr)
}
