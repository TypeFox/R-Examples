CreateTrueMean = function(tt,optns){
  # old mu_true
  
  tt[!(tt >= 0 & tt <= optns)] = 0
  mu =  (tt+sin(tt))
  return(mu[!is.na(mu)])
  
}
