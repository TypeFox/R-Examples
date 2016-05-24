## update mu and variance part
Step_2_Update_mu_var=function(model,num1)
{
  uu=0
  for(gg in (model$uniqueindex)){
    uu=uu+1
    if(gg <1){
      meanforupmu=(model$initial_mu_var$upvar[uu]*model$parameter$gamma0+model$parameter$eta0^2*sum(model$rstat[which(model$wholeindex==gg)]))/(model$initial_mu_var$upvar[uu]+model$parameter$eta0^2*sum(model$wholeindex==gg))
      varforupmu=model$initial_mu_var$upvar[uu]*model$parameter$eta0^2/(model$initial_mu_var$upvar[uu]+model$parameter$eta0^2*sum(model$wholeindex==gg))
      
      model$initial_mu_var$upmu[uu]<-rnorm(1,mean=meanforupmu,sd=sqrt(varforupmu))
      all0<-model$parameter$alpha0+sum(model$wholeindex==gg)*0.5
      laa0<-sum((model$rstat[which(model$wholeindex==gg)]-model$initial_mu_var$upmu[uu])^2)*0.5+model$parameter$beta0
      model$initial_mu_var$upvar[uu]<-rigamma(1,all0,laa0)
      
    }else{
      meanforupmu=(model$initial_mu_var$upvar[uu]*model$parameter$gamma1+model$parameter$eta1^2*sum(model$rstat[which(model$wholeindex==gg)]))/(model$initial_mu_var$upvar[uu]+model$parameter$eta1^2*sum(model$wholeindex==gg))
      varforupmu=model$initial_mu_var$upvar[uu]*model$parameter$eta1^2/(model$initial_mu_var$upvar[uu]+model$parameter$eta1^2*sum(model$wholeindex==gg))
      
      model$initial_mu_var$upmu[uu]<-rnorm(1,mean=meanforupmu,sd=sqrt(varforupmu))
      all1<-model$parameter$alpha1+sum(model$wholeindex==gg)*0.5
      laa1<-sum((model$rstat[which(model$wholeindex==gg)]-model$initial_mu_var$upmu[uu])^2)*0.5+model$parameter$beta1
      model$initial_mu_var$upvar[uu]<-rigamma(1,all1,laa1)
      
    }
    
  }
  
  if((model$sampleindex %in% model$originalindex)==0){
    if(model$sampleindex==min(model$myuniqueindex)){
      
      results=Gibbsfortheta(model$rstat[num1],model$parameter$beta0/(model$parameter$alpha0+0.5),0,model=model)
      newupmu=results[1]
      newupvar=results[2]
      model$initial_mu_var$upmu<-c(newupmu,model$initial_mu_var$upmu)
      model$initial_mu_var$upvar<-c(newupvar,model$initial_mu_var$upvar)
      
    }else if(model$sampleindex==max(model$myuniqueindex)){
      results=Gibbsfortheta(model$rstat[num1],model$parameter$beta1/(model$parameter$alpha1+0.5),1,model=model)
      newupmu=results[1]
      newupvar=results[2]
      model$initial_mu_var$upmu<-c(model$initial_mu_var$upmu,newupmu)
      model$initial_mu_var$upvar<-c(model$initial_mu_var$upvar,newupvar)
      
    }
  }
  return(model)
}
