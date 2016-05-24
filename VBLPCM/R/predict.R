predict.vblpcm<-function(object, ...)
  {
  d<-object$d
  N<-object$N
  XX_e<-object$XX_e
  V_xi_n<-object$V_xi_n
  V_xi_e<-object$V_xi_e
  V_psi2_n<-object$V_psi2_n
  V_psi2_e<-object$V_psi2_e
  V_z<-object$V_z
  V_sigma2<-object$V_sigma2
  DIST<-c(sqrt(as.matrix(dist(V_z, diag=1, upper=1)^2) + matrix(d*apply(expand.grid(V_sigma2,V_sigma2),1,sum),N))) 
  #DIST<-c(sqrt(as.matrix(-1/dist(V_z, diag=1, upper=1)^2) + matrix(d*apply(expand.grid(V_sigma2,V_sigma2),1,sum),N))) # NEW dists
  # the exponent of the expected log-likelihood is not the same as the expected likelihood
  cov1<-XX_e%*%V_xi_e-DIST
  cov2<-XX_e%*%V_psi2_e
  if (object$P_n > 0)
    {
    tmp2<-matrix(0,N,N)
    for (i in 1:N)
      for (j in 1:N)
        {
        if (object$model=="rsender")
          tmp2[i,j]=V_xi_n[i]
        if (object$model=="rreceiver")
          tmp2[i,j]=V_xi_n[j]
        if (object$model=="rsocial")
          tmp2[i,j]=V_xi_n[i,1]+V_xi_n[j,2]
	}
    cov1<-cov1+c(tmp2)
    tmp<-rep(V_psi2_n,N) 
    for (i in 1:N)
      for (j in 1:N)
        {
        if (object$model=="rsender")
          tmp2[i,j]=tmp[i]
        if (object$model=="rreceiver")
          tmp2[i,j]=tmp[j]
        if (object$model=="rsocial")
          tmp2[i,j]=tmp[i]+tmp[j]
	}
    cov2<-cov2+c(tmp2)
    }
  tmp<-cov1+0.5*cov2
  probs<-matrix(exp(tmp-log(1+exp(tmp))),N) 
  return(probs)
  }
