lmmlassoControl <-
function(tol=10^(-4),trace=1,maxIter=1000,maxArmijo=20,number=5,a_init=1,
                            delta=0.1,rho=0.001,gamma=0,lower=10^(-6),upper=10^8,seed=532,VarInt=c(0,10),CovInt=c(-5,5),thres=10^(-4))
  {

    list(tol=tol,trace=trace,maxIter=maxIter,maxArmijo=maxArmijo,number=number,
         a_init=a_init,delta=delta,rho=rho,gamma=gamma,lower=lower,upper=upper,seed=seed,VarInt=VarInt,CovInt=CovInt,thres=thres)
  }
