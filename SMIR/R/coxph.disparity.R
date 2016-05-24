## but need to account for smallest observations being censored!!!
coxph.disparity<-function(fit){
  if (class(fit)!="coxph") stop("Argument needs to be class ``coxph''")
  require(survival)
bhz<-basehaz(fit,centered=TRUE)
bhz$haz<-diff(c(0,bhz$hazard))
bhz$tint<-diff(c(0,bhz$time))
ti <- apply(outer(fit$y[,1],bhz$time,">="),1,sum)
  -2*(sum(fit$y[,2]*(ti!=0)*(log(bhz$haz[ti+(ti==0)]/bhz$tint[ti+(ti==0)])+fit$linear)-(ti!=0)*bhz$hazard[ti+(ti==0)]*exp(fit$linear)))
}
