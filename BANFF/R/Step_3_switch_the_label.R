## switch the label
Step_3_switch_the_label=function(model,num1)
{
  model$uniqueindex<-sort(unique(model$wholeindex))
  test2whole=model$wholeindex
  test2mu=model$initial_mu_var$upmu
  test2var=model$initial_mu_var$upvar
  
  if(length(model$uniqueindex)>1){
    orderindex=rank(model$initial_mu_var$upmu)
    ssupmu=0
    ssupvar=0
    sswholeindex<-rep(NA,length(model$rstat))
    for(kk in 1:length(model$uniqueindex)){
      ssupmu[kk]=model$initial_mu_var$upmu[which(orderindex==kk)]
      ssupvar[kk]=model$initial_mu_var$upvar[which(orderindex==kk)]
      dudu=which(orderindex==kk)
      sswholeindex[which(model$wholeindex==model$uniqueindex[dudu])]=model$uniqueindex[kk]
      
      
    }
    model$initial_mu_var$upmu=ssupmu
    model$initial_mu_var$upvar=ssupvar
    model$wholeindex=sswholeindex
  }
  return(model)
}
