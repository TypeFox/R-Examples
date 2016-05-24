extractplate = function(datbefore, datafter, plate, replicate){
  datbefore = datbefore[[replicate]]
  datafter = datafter[[replicate]]

  if (plate == 1){
    datbefore = datbefore[seq(1,nrow(datbefore),2),]
    datbefore = datbefore[,seq(1,24,2)]

    datafter = datafter[seq(1,nrow(datafter),2),]
    datafter = datafter[,seq(1,24,2)]
  }
  else if (plate == 2){
    datbefore = datbefore[seq(1,nrow(datbefore),2),]
    datbefore = datbefore[,seq(2,24,2)]

    datafter = datafter[seq(1,nrow(datafter),2),]
    datafter = datafter[,seq(2,24,2)]
  }
  else if (plate == 3){
    datbefore = datbefore[seq(2,nrow(datbefore),2),]
    datbefore = datbefore[,seq(1,24,2)]

    datafter = datafter[seq(2,nrow(datafter),2),]
    datafter = datafter[,seq(1,24,2)]
  }
  else if (plate == 4){
    datbefore = datbefore[seq(2,nrow(datbefore),2),]
    datbefore = datbefore[,seq(2,24,2)]

    datafter = datafter[seq(2,nrow(datafter),2),]
    datafter = datafter[,seq(2,24,2)]
  }
  else 
    stop ("unknown plate.")

  datall = list(datbefore=datbefore, datafter=datafter)
  return(datall)
}
