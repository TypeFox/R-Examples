mix.qn <-
function(y0,ctrl.id,NQ,mix.param,max.q=0.95,low){ #anything greater than max.q, use param
  ##
  ##NQ:normalized quantiles
  ##
  ECDF=ecdf(y0[ctrl.id])
  Q=ECDF(y0)
  id0=which(Q<0.05) ##lower tail, set to lower bound
  id2=which(Q>max.q) ##higher tail, use parametric version
  B=length(id2)
  Q2=max.q+(1-max.q)*(1:B)/(B+1)
    
  y2o=order(y0[id2])
  Q[id2][y2o]=Q2
  ## qs[id0]=0.05
  
  ## < 0.05, set to 0.05; greater than .95, use normal.mixture;between, use empirical
  ynorm=vector(mode="numeric",length=length(y0))
  ynorm[id0]=low
  ynorm[-c(id0,id2)]=quantile(NQ,Q[-c(id0,id2)])
  ynorm[id2]=qnorMix(Q[id2],mix.param)
  ynorm
}

