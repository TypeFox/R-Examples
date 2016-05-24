######DMP1 Iteration functions
Step1_1_Update_gi_zi=function(model,num1)
{
  model$uniqueindex<-sort(unique(model$wholeindex))
  model$originalindex<-model$uniqueindex
  sampleprob<-rep(NA,length(model$uniqueindex)+2)
  
  
  idx = which(model$net[num1,]==1)
  
  pro0=sum((model$wholeindex[idx]<1))
  pro1=sum((model$wholeindex[idx]>0))
  ## calculate mk
  propro0=sum(model$wholeindex<1)
  propro1=sum(model$wholeindex>0)
  
  
  sumforg=0
  uu<-0
  for(gg in model$uniqueindex){
    uu<-uu+1
    sumforg[uu]<-sum(model$wholeindex==gg)-(model$wholeindex[num1]==gg)
    if(gg<1){
      sampleprob[uu+1]=log(model$parameter$pi0)+(2*model$parameter$rho0*pro0-(model$rstat[num1]-model$initial_mu_var$upmu[uu])^2*0.5/model$initial_mu_var$upvar[uu])-0.5*log(model$initial_mu_var$upvar[uu])+log(sumforg[uu])-log(model$parameter$tau0+propro0-1)
      
    }else{
      sampleprob[uu+1]=log(1-model$parameter$pi0)+(2*model$parameter$rho1*pro1-(model$rstat[num1]-model$initial_mu_var$upmu[uu])^2*0.5/model$initial_mu_var$upvar[uu])-0.5*log(model$initial_mu_var$upvar[uu])+log(sumforg[uu])-log(model$parameter$tau1+propro1-1)
    }
  }
  
  test1whole=model$wholeindex
  
  sampleprob[1]= log(model$parameter$pi0)+(2*model$parameter$rho0*pro0)+log(model$parameter$tau0)-log(model$parameter$tau0+propro0-1)+log(Innerfunc(0,model$rstat[num1],model=model))
  sampleprob[length(model$uniqueindex)+2]= log(1-model$parameter$pi0)+(2*model$parameter$rho1*pro1)+log(model$parameter$tau1)-log(model$parameter$tau1+propro1-1)+log(Innerfunc(1,model$rstat[num1],model=model))
  newsampleprob<-0
  for(uu in 1:length(sampleprob)){ newsampleprob[uu]=1/sum(exp(sampleprob-sampleprob[uu]))}
  newsampleprob[which(newsampleprob=="NaN")]=0
  model$origindex<-model$wholeindex[num1]
  model$wholeindex[num1]=sample((min(model$wholeindex)-1):(max(model$wholeindex)+1),1,prob=newsampleprob)
  model$sampleindex=model$wholeindex[num1]  
  return(model)
}
