pba <-
function(x, tr, newx, k=1, global=0.5){

  if(!is.data.frame(x))             stop('x must be a data.frame')
  if(length(tr) != nrow(x))         stop('tr length not equal to x rows')
  if(!is.data.frame(newx))          stop('newx must be a data.frame')
  if(ncol(newx) != ncol(x))         stop('newx does not have same 
                                         number of covariates as x')
  if(prod(names(newx)!= names(x)))  stop('newx covariates differ from x')
  if(nrow(newx) > 1)                stop('newx has more than one row')
  if(k < 0)                         stop('k must be 0, Inf, 
                                         or a ratio of positive odd integers')
  if(global <= 0 | global >= 1)     stop('global is proportion to be treated: 
                                         need 0 < global < 1')

  ## obtain fitted propensity score
  tmodel = glm(tr ~., family=binomial('logit'), data=x)
  phat = predict(tmodel, newdata=newx, type='response')[[1]]

  ## deterministic PBA
  if(k == Inf){
    if(phat < global){
      ptreat = 1
      newtr = 1
    }else if(phat == global){
      ptreat = global
      newtr = rbinom(1, size=1, prob=ptreat)
    }
    else if(phat > global){
      ptreat = 0
      newtr = 0
    }
  ## PBA with k = 0 (approx. pure randomization)
  }else if(k == 0){
    if(phat == 0){
      ptreat = 1
      newtr = 1
    }else if(phat == 1){
      ptreat = 0
      newtr = 0
    }
    else{
      ptreat = 0.5
      newtr = rbinom(1, size=1, prob=ptreat)
    }
  ## random PBA
  }else{
    ptreat = piFunction(fit=phat, kparam=k, qparam=global)
    newtr = rbinom(1, size=1, prob=ptreat)
  }

  inp = list(x=x, tr=tr, newx=newx, k=k, global=global)
  res = list(phat = phat, ptreat=ptreat, newtr=newtr)
  return(list(results=res, input=inp))
}
