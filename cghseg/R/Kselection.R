Kselection <- function(Y,t.est,loglik,Kseq,CGHo){

  Kmax  = max(Kseq)
  J     = -loglik
  Kh    = Kmax
  n     = length(Y)
  SSall = sum((Y-mean(Y,na.rm=TRUE))^2,na.rm=TRUE)
  BIC   = rep(0,length(Kseq))
  n     = length(Y)

  for (k in Kseq){
    th      = t.est[k,1:k]
    rupt    = matrix(ncol=2,c(c(1,th[1:k-1]+1),th))
    mu      = data.frame(begin = rupt[,1],
      end   = rupt[,2],
      mean  = apply(rupt,1,FUN=function(z) mean(Y[z[1]:z[2]], na.rm=TRUE)))
    gamma.coef = lgamma(0.5*(n-k-1)) - lgamma(0.5*(n+1))
    nk   = mu$end-mu$begin+1
    SSbg = sum(nk*(mu$mean-mean(Y,na.rm=TRUE))^2,na.rm=TRUE)
    if (SSbg ==0){
      BIC[k] =  gamma.coef + 0.5*k*log(SSall)-0.5 * sum(log(nk))+(0.5-k)*log(n)           
    } else {
      BIC[k] = 0.5*(n-k+1) * log(1+ (SSbg)/(SSall-SSbg)) + gamma.coef + 0.5*k*log(SSall)-0.5 * sum(log(nk))+(0.5-k)*log(n)                       }
  } # end k
  Kh         = Kseq[which.max(BIC)]
  invisible(Kh)
}#end fonction




