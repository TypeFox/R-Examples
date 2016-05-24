"getStdErrClp" <- function(group,  multimodel, thetalist,
                           clpindepX, rlist, rawtheta) {

  QR.temp <- rlist$QR 
  m <- multimodel@modellist
  dset <- group[[1]][2]
  clpind <- group[[1]][1]
  cp <- multimodel@fit@resultlist[[dset]]@cp[[clpind]]
  sigma_2 <- multimodel@fit@nlsres$sumonls$sigma^2
  if(class(multimodel@fit@nlsres$onls) == "timp.nls") {
    R <- multimodel@fit@nlsres$onls$m$Rmat()
    R_inv <- chol2inv(R) 
  }
  else{ 
    R <- multimodel@fit@nlsres$onls$hessian
    R_inv <- solve(R)
  }
  R_X <- qr.R(QR.temp)
  R_X_inv <- backsolve(R_X,diag(ncol(R_X)))
  X_pseudo <- tcrossprod(R_X_inv, qr.Q(QR.temp))
  numenv <- new.env()
  assign("model", m[[dset]], envir = numenv)
  assign("group", group, envir = numenv)
  assign("rawtheta", rawtheta, envir = numenv)
  assign("multimodel", multimodel, envir = numenv)
  assign("thetalist", thetalist, envir = numenv)
  assign("clpindepX", clpindepX, envir = numenv)
  assign("finished", FALSE, envir = numenv)
  assign("returnX", TRUE, envir = numenv)
  if(m[[ dset ]]@clpdep) 
    s_e <- body(selectMethod("residPart", m[[ dset ]]@mod_type)@.Data)
  else s_e <- body(selectMethod("residPart", m[[ dset ]]@mod_type)@.Data)
  X_gradient <- attr( numericDeriv(expr=s_e, rho = numenv, 
                                   theta = c("rawtheta")), "gradient")
  dim(X_gradient) <- c(dim(QR.temp$qr),length(rawtheta)) 
  G <- matrix(nrow=length(cp), ncol = length(rawtheta)) 
  for(i in 1:length(rawtheta)) 
    G[,i] <- X_pseudo %*%  ( as.matrix(X_gradient[,,i]) %*% cp) 
  G_R_inv <- G %*% R_inv
  ## if A R = G, R^T A^T = G^T 
  ## G R_inv R_inv^T G^T 
  ## A A^T  
  A_T <- solve(t(R), t(G))
  Bloc1 <- tcrossprod(X_pseudo, X_pseudo) 
  std_err_clp <- sqrt( sigma_2 * diag(Bloc1 + crossprod(A_T, A_T)))
  for (i in 1:length(group)) 
    multimodel@fit@resultlist[[group[[i]][2]]]@std_err_clp[[group[[i]][1]]] <- std_err_clp 
  multimodel 
}
