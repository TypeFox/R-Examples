surv.exp.gt.model <- function(pilm,lm,gtprev,GRR,zmodel,interval)
  {
    ww=GRR^zmodel
    foo=function(lam0) {
        ww=-lam0*ww*lm
        pilm-sum(exp(ww)%*%gtprev)
    }
    ww*uniroot(foo,interval)$root
  }
