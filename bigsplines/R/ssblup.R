ssblup <-
  function(Zty,ZtX,ZtZ,coefs,rdm,tau){
    
    Zte <- Zty-ZtX%*%coefs
    if(length(tau)==1L){
      if(is.matrix(ZtZ)){
        bhat <- as.numeric(tau*(diag(nrow(ZtZ))-ZtZ%*%pinvsm(diag(rep(1/tau,nrow(ZtZ)))+ZtZ))%*%Zte)
      } else {
        bhat <- as.numeric(tau*(1-ZtZ/((1/tau)+ZtZ))*Zte)
      }
    } else {
      if(is.matrix(ZtZ)){
        bhat <- as.numeric(rep(tau,rdm)*(diag(nrow(ZtZ))-ZtZ%*%pinvsm(diag(rep(1/tau,rdm))+ZtZ))%*%Zte)
      } else {
        bhat <- as.numeric(rep(tau,rdm)*(1-ZtZ/((rep(1/tau,rdm))+ZtZ))*Zte)
      }
    }
    names(bhat) <- names(Zty)
    return(bhat)
    
  }