PKselection <- function(Y,Kseq,param,loglik,CGHo){
  
  Kmin         = min(Kseq)
  Kmax         = max(Kseq)
  P            = CGHo["nblevels"]
  S            = 0.5
  present.data = which(!is.na(Y))
  missing.data = which(is.na(Y))
  x            = Y[present.data]
  lmax         = floor(CGHo["lmax"]*length(x))
  n            = length(x)

  ##################################################
  #                                                #
  #        MODIFIED BIC FOR MODEL SELECTION        #
  #                                                #
  ##################################################
  
  SSall      = sum((x-mean(x))^2)
  mBIC       = rep(-Inf,Kmax)  
  gamma.coef = lgamma(0.5*(n-P-1)) - lgamma(0.5*(n-1))

  
  if (P==1){
    mBIC = gamma.coef+0.5*log(SSall) - Kseq*log(n)
  } else {
    for (k in (P:Kmax)){          
      SSbg       = 0
      np         = c()
      phitmp     = param[[k]]$phi
      rupt       = param[[k]]$rupt
      
      nk         = rupt[,2]-rupt[,1]+1
      logdensity = t( apply(rupt,1,FUN=function(y) logdens(   x[ y[1]:y[2] ] ,P,phitmp)))          
      tau        = Estep(logdensity,phitmp)
      pop        = apply(tau,1,which.max)
      cluster    = c()      
      for ( j in (1:k) ) {
        cluster[rupt[j,1]:rupt[j,2]] = rep(pop[j],rupt[j,2]-rupt[j,1]+1)
      }

      for (ell in c(1:P)){
        np[ell] = sum(nk[pop==ell])
        if (np[ell] == 0){
          np[ell] = 1
        }
        mp      = mean(x[cluster==ell])
        if (is.na(mp)){
          mp=0
        }
        SSbg    = SSbg+np[ell]*(mp-mean(x))^2        
      }
      mBIC[k] = (n-P-1)/2 * log(1+ (SSbg)/(SSall-SSbg)) + gamma.coef + 
        0.5*P*log(SSall)-0.5 * sum(log(np))-(k-0.5)*log(n)                       
    } # end k
    mBIC[mBIC==0]     = -Inf
    mBIC[is.na(mBIC)] = -Inf
  }
  Kselect           = which.max(mBIC)
  invisible(Kselect)
  
}
