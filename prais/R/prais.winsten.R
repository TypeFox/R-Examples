prais.winsten <-
  function(formula,data,iter=50,rho=0,tol=.00000001){
    mod<-model.frame(formula,data=data)
    lm<-lm(mod)
    n<-length(mod[,1])
    list.rho<-c(0)
    imax<-ncol(mod)-1
    fo<-as.formula(paste("y ~ -1 + x0 +",paste(paste0("x",1:imax), collapse= "+")))
    # In case rho is known
    if (rho!=0) {
      y<-c((1-rho^2)^(1/2)*mod[1,1],mod[2:n,1]-rho*mod[1:(n-1),1])
      x0<-c((1-rho^2)^(1/2),rep(1-rho,n-1))
      for (i in 1:imax) {
        x<-c((1-rho^2)^(1/2)*mod[1,(i+1)],mod[2:n,(i+1)]-rho*mod[1:(n-1),(i+1)])
        assign(paste("x",i,sep=""),x)
      }
      lm<-lm(fo)
      j<-1
      rho.tstat<-"none"
    }
    # Include the estimation of rho (usual case)
    else {
      res<-lm$res
      res_1<-c(NA,res[-n])
      for (i in 1:iter){
        rho.lm<-lm(res ~ res_1 -1)
        rho<-rho.lm$coeff[1]
        if (abs(rho-tail(list.rho,n=1))<tol) {
          j<-i
          break
        }
        else {
          list.rho<-append(list.rho,rho)
          y<-c((1-rho^2)^(1/2)*mod[1,1],mod[2:n,1]-rho*mod[1:(n-1),1])
          x0<-c((1-rho^2)^(1/2),rep(1-rho,n-1))
          for (k in 1:imax) {
            x<-c((1-rho^2)^(1/2)*mod[1,(k+1)],mod[2:n,(k+1)]-rho*mod[1:(n-1),(k+1)])
            assign(paste("x",k,sep=""),x)
          }
          lm<-lm(fo)
          fit<-as.vector(rep(lm$coef[1],n))+as.vector(as.matrix(mod[,2:(imax+1)]) %*% lm$coef[2:(imax+1)])
          res<-mod[,1]-fit
          res_1<-c(NA,res[-n])
          j<-i
          rho.tstat<-summary(rho.lm)$coef[1,3]
        }
      }
    }
    if (iter==50) i<-j-1 else i<-j
    attr(lm$coefficients,"names")<-c('Intercept',names(mod)[2:ncol(mod)])
    s<-summary(lm)
    r<-data.frame('Rho'=rho,'Rho.t-statistic'=rho.tstat,'Iterations'=i,row.names=c(''))
    results<-list(s,r)
    return(results)
  }