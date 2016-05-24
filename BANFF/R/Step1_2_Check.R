## check whether one cluster is disappear
Step1_2_Check=function(model,num1)
{
  model$myuniqueindex<-sort(unique(model$wholeindex))
  if(sum(model$wholeindex==model$origindex)==0){
    
    deletindex=which(model$originalindex==model$origindex)
    
    model$initial_mu_var$upmu=model$initial_mu_var$upmu[-deletindex]
    model$initial_mu_var$upvar=model$initial_mu_var$upvar[-deletindex]
    if(model$origindex<1){
      model$wholeindex[which(model$wholeindex<model$origindex)]=model$wholeindex[which(model$wholeindex<model$origindex)]+1
    }else{
      
      model$wholeindex[which(model$wholeindex>model$origindex)]=model$wholeindex[which(model$wholeindex>model$origindex)]-1
      
    }
    model$uniqueindex<-sort(unique(model$wholeindex))
    if((model$sampleindex %in% model$originalindex)==0){
      if(model$sampleindex<1){model$uniqueindex=model$uniqueindex[-1]}else{model$uniqueindex=model$uniqueindex[-length(model$uniqueindex)]}
      
    }
    
    
  }
  
  
  test1mu=model$initial_mu_var$upmu
  test1var=model$initial_mu_var$upvar
  test1unique=model$uniqueindex
  return(model)
}
