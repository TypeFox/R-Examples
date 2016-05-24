CalculScores <-function(sequencia, logodds){
  score <-0

  for(i in 1:dim(logodds)[1]){
    if(sequencia[i] == "A"){
      score <- logodds[i,1]+score
    }else if(sequencia[i]== "T"){
      score <- logodds[i,2]+score
    }else if(sequencia[i]== "C"){
      score <- logodds[i,3]+score
     }else if(sequencia[i]== "G"){
      score <- logodds[i,4]+score
     }
     
  }
   score
}

