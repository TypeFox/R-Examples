vblpcmKL<-function(x)
  {
  model<-x$model
  P_n<-x$P_n
  P_e<-x$P_e
  d<-x$d
  N<-x$N
  NE<-x$NE
  NnonE<-x$NnonE
  NM<-x$NM
  G<-x$G
  Y<-x$Y
  E<-x$E
  nonE<-x$nonE
  M<-x$M
  numedges<-x$numedges
  EnonE<-x$EnonE
  diam<-x$diam
  hopslist<-x$hopslist
  XX_e<-x$XX_e
  V_xi_n<-x$V_xi_n
  V_xi_e<-x$V_xi_e
  V_psi2_n<-x$V_psi2_n
  V_psi2_e<-x$V_psi2_e
  V_z<-x$V_z
  V_sigma2<-x$V_sigma2
  V_eta<-x$V_eta
  V_lambda<-x$V_lambda
  V_omega2<-x$V_omega2
  V_nu<-x$V_nu
  V_alpha<-x$V_alpha
  xi<-x$xi
  psi2<-x$psi2
  sigma02<-x$sigma02
  omega2<-x$omega2
  nu<-x$nu
  alpha<-x$alpha
  inv_sigma02<-x$inv_sigma02
  seed<-x$seed
  NC=x$NC
  KL=0
  imodel=switch(model, plain=0, rsender=1, rreceiver=2, rsocial=3)
  total_KL<-function(imodel, P_n, P_e, D, N, NE, NnonE, NM, G, Y, E, nonE, M, numedges, EnonE, diam, hopslist, XX_e, 
                     V_xi_n, V_xi_e, V_psi2_n, V_psi2_e, V_z, V_sigma2, V_eta, V_lambda, V_omega2, V_nu, V_alpha, xi, 
                     psi2, sigma02, omega2, nu, alpha, inv_sigma02, seed, NC, KL) 
  		     {
                     ans<-.C("KL_total", NAOK=TRUE, 
                     imodel=as.integer(imodel), P_n=as.integer(P_n),P_e=as.integer(P_e),D=as.integer(d), N=as.integer(N), 
  		     NE=as.integer(NE), NnonE=as.integer(NnonE), NM=as.integer(NM),
                     G=as.integer(G), Y=as.numeric(t(Y)), E=as.integer(t(E)), nonE=as.integer(t(nonE)), M=as.integer(t(M)),
                     numedges=as.integer(t(numedges)), EnonE=as.integer(t(EnonE)),
		     diam=as.integer(diam), hopslist=as.integer(t(hopslist)), 
                     XX_e=as.double(t(XX_e)), V_xi_n=as.double(V_xi_n), V_xi_e=as.double(V_xi_e), 
                     V_psi2_n=as.double(V_psi2_n), V_psi2_e=as.double(V_psi2_e), V_z=as.double(t(V_z)),
                     V_sigma2=as.double(V_sigma2), V_eta=as.double(t(V_eta)), V_lambda=as.double(t(V_lambda)),
                     V_omega2=as.double(V_omega2), V_nu=as.double(V_nu), V_alpha=as.double(V_alpha),
                     xi=as.double(xi), psi2=as.double(psi2), sigma02=as.double(sigma02),
                     omega2=as.double(omega2), nu=as.double(nu), alpha=as.double(alpha),
                     inv_sigma02=as.double(inv_sigma02), seed=as.double(seed),
  		     NC=as.integer(NC), KL=as.double(KL), PACKAGE="VBLPCM")
                     return(ans)
                     }
  final_KL<-total_KL(imodel, P_n, P_e, d, N, NE, NnonE, NM, G, Y, E, nonE, M, numedges, EnonE, diam, hopslist, XX_e, V_xi_n, V_xi_e, 
                     V_psi2_n, V_psi2_e, V_z, V_sigma2, V_eta, V_lambda, V_omega2, V_nu, V_alpha, xi, psi2, sigma02, omega2, nu, alpha,
  		     inv_sigma02, seed, NC, KL)
  cat("KL distance to true posterior is ", final_KL$KL, "+ constant \n")
  final_KL$KL
  }

