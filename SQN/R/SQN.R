SQN <-
function(y,N.mix=5,ctrl.id,model.weight=.9){
  QE=apply(y[ctrl.id,],2,sort)
  QN=apply(QE,1,median)
  mix.param=Mclust(QN,G=N.mix)$parameters
  mix.param=norMix(mu=mix.param$mean,sig2=mix.param$variance$sigmasq,w=mix.param$pro)

  qq=seq(1/(2*length(QN)),1-1/(2*length(QN)),1/length(QN))
  qq=qnorMix(qq,mix.param)
  QN1=QN*(1-model.weight)+qq*model.weight
  
  ynorm=apply(y,2,mix.qn,ctrl.id,NQ=QN1,mix.param=mix.param,max.q=0.95,low=quantile(QN1,.05))
  
}

