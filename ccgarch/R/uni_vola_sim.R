# simulating univariate volatility
   uni.vola.sim <- function(a, nobs, d.f=Inf, cut=1000){
      nobs <- nobs + cut
#      cond.h <- .Call("uni_vola_sim", nobs, a, d.f, PACKAGE="ccgarch")
      cond.h <- .Call("uni_vola_sim", nobs, a, d.f)
      list( h = cond.h[[1]][(cut+1):nobs],
            eps = cond.h[[2]][(cut+1):nobs])
   }
