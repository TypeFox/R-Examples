fitness_count <-
function(Al,para,Data){

  ### as from our exploration on genetic algorithm, we should have more population, perform crossover, just forget about mutation ###
  R<-aligned_string(Data,para)
  #R<-alignment(Data,para)
  
  fitness<-vector()
  
  for(i in 1:length(para[,1])){
    fitness[i]<-sum(R[[1]][i,]==Al[1,]&R[[2]][i,]==Al[2,])
  }
  
  return(fitness)
  
}
