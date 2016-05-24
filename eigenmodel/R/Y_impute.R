"Y_impute" <-
function() {

  sd_zq<-1/sqrt(pp_zq)
  zq<-c(-Inf,rep(NA,max(Ranks,na.rm=TRUE)-1),Inf)
  for(ry in 1:(max(Ranks,na.rm=TRUE)-1)){
    ub<-suppressWarnings(min(Z[ Ranks==ry+1 ],na.rm=TRUE ) )
    lb<-suppressWarnings(max(Z[ Ranks==ry ],na.rm=TRUE ) )
    zq[ry+1]<-  qnorm( runif(1,pnorm(lb,0,sd_zq),pnorm(ub,0,sd_zq)),0,sd_zq  )
                                     }

  zhat<- Z[upper.tri(Z) & is.na(Y)]
  lb<-outer(zhat,zq[-1],"<")
  ub<-outer(zhat,zq[-length(zq)],">")
  est<-lb & ub
  yhat<- est%*% sort(unique(c(Y)))

  Y[ upper.tri(Z) & is.na(Y) ] <- yhat  
  Y[ lower.tri(Y) ]<- 0 
  Y+t(Y)             }

