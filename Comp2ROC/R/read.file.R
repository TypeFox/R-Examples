read.file <-
function(name.file.csv,header.status=TRUE,separator=";",decimal=",", modality1, testdirection1, modality2,testdirection2,status1,related=TRUE,status2=NULL){
    dataFromFile = read.csv(name.file.csv,header=header.status,sep=separator,dec=decimal)    
    if (is.null(status2)) {
      dat = dataFromFile 
    } else {
      dat = list()
      dat[[modality1]] = dataFromFile[[modality1]]
      dat[[modality2]] = dataFromFile[[modality2]]
      dat[[status1]] = dataFromFile[[status1]]
      dat[[status2]] = dataFromFile[[status2]]
    }
    
    #Select modality
    if(testdirection1){
      mod1.ind=dat[modality1]
    } else{
      mod1.ind=(-1)*dat[modality1]
    }
    if(testdirection2){
      mod2.ind=dat[modality2]
    } else {
      mod2.ind=(-1)*dat[modality2]
    }
    
    sim1.ind=mod1.ind
    sim2.ind=mod2.ind
    
    if(related){
      sim1.sta=dat[status1]
      sim2.sta=sim1.sta
    } else{
      sim1.sta=dat[status1]
      sim2.sta=dat[status2]
    }
    
    #Clean NAs
    if (!is.null(status2)) {
      for (i in 1:length(sim1.ind[[modality1]])){
        if (is.na(sim1.ind[[modality1]][i])) {
          indexToRemove = i-1
          sim1.ind[[modality1]] = sim1.ind[[modality1]][0:indexToRemove]
          break
        }
      }
      for (i in 1:length(sim1.sta[[status1]])){
        if (is.na(sim1.sta[[status1]][i])) {
          indexToRemove = i-1
          sim1.sta[[status1]] = sim1.sta[[status1]][0:indexToRemove]
          break
        }
      }
      for (i in 1:length(sim2.ind[[modality2]])){
        if (is.na(sim2.ind[[modality2]][i])) {
          indexToRemove = i-1
          sim2.ind[[modality2]] = sim2.ind[[modality2]][0:indexToRemove]
          break
        }
      }
      for (i in 1:length(sim2.sta[[status2]])){
        if (is.na(sim2.sta[[status2]][i])) {
          indexToRemove = i-1
          sim2.sta[[status2]] = sim2.sta[[status2]][0:indexToRemove]
          break
        }
      }
    }
    
    answer=list(sim1.ind,sim2.ind,sim1.sta,sim2.sta)
    
    return(answer)
  }
