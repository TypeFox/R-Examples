
fit_fun = function(ImpCov,SampCov,Areg,lambda,alpha,type,pen_vec){


  m = dim(ImpCov)[1]
  IntCol = which(colnames(Areg) == "1")

  if(type=="none"){

    fit = log(det(ImpCov)) + trace(SampCov %*% solve(ImpCov)) - log(det(SampCov))  - m
  }

  else if(type=="ridge"){
   #ridge
  fit = log(det(ImpCov)) + trace(SampCov %*% solve(ImpCov)) - log(det(SampCov)) - m  + 2*lambda * sum(pen_vec)
  }

  else if(type=="lasso"){
  #lasso
  fit = log(det(ImpCov)) + trace(SampCov %*% solve(ImpCov)) - log(det(SampCov)) - m  + 2*lambda * sum(abs(pen_vec))
  }

  else if(type=="diff"){
    #lasso
    fit = log(det(ImpCov)) + trace(SampCov %*% solve(ImpCov)) - log(det(SampCov)) - m  + 2*lambda * sd(Areg[Areg != 0])
  }
  else if(type=="diff_growth"){
    #lasso
    fit = log(det(ImpCov)) + trace(SampCov %*% solve(ImpCov)) - log(det(SampCov)) - m +
          2*lambda * sum(abs(Areg[,ncol(Areg)] - c(0,1,2,3,0,0,0)))
  }
 # else if(type=="gLasso"){
    #lasso
  #  fit = log(det(solve(ImpCov))) + trace(SampCov %*% solve(solve(ImpCov)))  - lambda * sum(abs(solve(ImpCov)))
 # }
  else if(type=="gLasso"){
    #lasso
    imp_pen <- solve(ImpCov)
    diag(imp_pen) <- 0
    imp_pen2 <- sum(abs(imp_pen))
    fit = log(det(solve(ImpCov))) + trace(SampCov %*% solve(ImpCov))  + 2*lambda * imp_pen2
  }
  else if(type=="gRidge"){
    imp_pen <- solve(ImpCov)
    diag(imp_pen) <- 0
    imp_pen2 <- sum(imp_pen %*% t(imp_pen))
    fit = log(det(ImpCov)) + trace(SampCov %*% solve(ImpCov))  + 2*lambda * imp_pen2
  }
  else if(type=="enet"){
    #elastic net
    fit = log(det(ImpCov)) + trace(SampCov %*% solve(ImpCov)) - log(det(SampCov))  +
      2*lambda * sum(alpha*(Areg*Areg)  + (1- alpha)*abs(Areg))
  }
  else if(type=="ols_lasso"){

    fit = 0.5 * trace((SampCov - ImpCov)^2)  + 2*lambda * sum(abs(Areg))
  }

  0.5*fit


  # ------------------------- penalized fit ------------------------------
  #lambda = -0.15

  #fit.reg = log(det(ImpCov)) + trace(cov %*% solve(ImpCov)) - log(det(cov)) - m + lambda*cov

  #fit.reg


}
