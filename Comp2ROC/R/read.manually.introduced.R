read.manually.introduced <-
function(dat,modality1,testdirection1,modality2,testdirection2,status1,related=TRUE,status2=NULL){
  #Select modality
  if(testdirection1){
    mod1.ind=dat[modality1]
  } else{
    mod1.ind=(-1)*as.numeric(unlist(dat[modality1]))
  }
  if(testdirection2){
    mod2.ind=dat[modality2]
  } else{
    mod2.ind=(-1)*as.numeric(unlist(dat[modality2]))
  }
  
  sim1.ind=mod1.ind
  sim2.ind=mod2.ind
  
  if(related){
    sim1.sta=dat[status1]
    sim2.sta=sim1.sta
  }
  else{
    sim1.sta=dat[status1]
    sim2.sta=dat[status2]
  }
  answer=list(sim1.ind,sim2.ind,sim1.sta,sim2.sta)
  return(answer)
  
}
