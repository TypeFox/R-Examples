vblpcmfit<-function(variational.start, STEPS=50, maxiter=100, tol=1e-6, NC=NULL, seed=NaN, d_vector=rep(TRUE,9))
  {
  if (length(d_vector)!=9)
    stop("You must supply a d_vector of length 9. Please refer to the help file for vblpcmfit\n")
  if (!is.null(NC))
    {
    NC<-as.integer(NC)
    if (NC==0)
      stop("Cannot use zero controls per case\n")
    cat("Using ", NC, "controls per case in case-control sampler\n")
    }
  if (is.nan(seed))
    seed=runif(1,0,1e6) # set the seed
  set.seed(seed) # use this to seed the random number generator in R
  P_n<-variational.start$P_n
  P_e<-variational.start$P_e
  model<-variational.start$model
  d<-variational.start$d
  N<-variational.start$N
  NE<-variational.start$NE
  NnonE<-variational.start$NnonE
  if (is.null(NC))
    NC<-NnonE # use all non-edges
  NM<-variational.start$NM
  G<-variational.start$G
  Y<-variational.start$Y
  E<-variational.start$E
  nonE<-variational.start$nonE
  M<-variational.start$M
  numedges<-variational.start$numedges
  EnonE<-variational.start$EnonE
  diam<-variational.start$diam
  hopslist<-variational.start$hopslist
  XX_e<-variational.start$XX_e
  V_xi_n<-variational.start$V_xi_n
  V_xi_e<-variational.start$V_xi_e
  V_psi2_n<-variational.start$V_psi2_n
  V_psi2_e<-variational.start$V_psi2_e
  V_z<-variational.start$V_z
  V_sigma2<-variational.start$V_sigma2
  V_eta<-variational.start$V_eta
  V_lambda<-variational.start$V_lambda
  V_omega2<-variational.start$V_omega2
  V_nu<-variational.start$V_nu
  V_alpha<-variational.start$V_alpha
  xi<-variational.start$xi
  psi2<-variational.start$psi2
  sigma02<-variational.start$sigma02
  omega2<-variational.start$omega2
  nu<-variational.start$nu
  alpha<-variational.start$alpha
  inv_sigma02<-variational.start$inv_sigma02
  imodel=switch(model, plain=0, rsender=1, rreceiver=2, rsocial=3)
  conv=0 # not converged to start with
  out<-.C("Rf_VB_bbs", NAOK=TRUE, imodel=as.integer(imodel), steps=as.integer(STEPS), max_iter=as.integer(maxiter), P_n=as.integer(P_n), 
           P_e=as.integer(P_e), D=as.integer(d), N=as.integer(N), NE=as.integer(NE), NnonE=as.integer(NnonE), NM=as.integer(NM),
           G=as.integer(G), Y=as.numeric(t(Y)), E=as.integer(t(E)), nonE=as.integer(t(nonE)), M=as.integer(t(M)),
  	   numedges=as.integer(t(numedges)), EnonE=as.integer(t(EnonE)), diam=as.integer(diam),
	   hopslist=as.integer(t(hopslist)), XX_e=as.double(t(XX_e)),
           V_xi_n=as.double((V_xi_n)), V_xi_e=as.double(V_xi_e), V_psi2_n=as.double(V_psi2_n), 
           V_psi2_e=as.double(V_psi2_e), V_z=as.double(t(V_z)), V_sigma2=as.double(V_sigma2), 
           V_eta=as.double(t(V_eta)), V_lambda=as.double(t(V_lambda)),
           V_omega2=as.double(V_omega2), V_nu=as.double(V_nu), V_alpha=as.double(V_alpha),
           xi=as.double(xi), psi2=as.double(psi2), sigma02=as.double(sigma02),
           omega2=as.double(omega2), nu=as.double(nu), alpha=as.double(alpha),
           inv_sigma02=as.double(inv_sigma02), tol=as.double(tol), NC=as.integer(NC), 
	   seed=as.double(seed), d_vector=as.double(d_vector), conv=as.integer(conv),
	   PACKAGE="VBLPCM")
  if (model=="plain") 
    V_xi_n<-NaN
  if (model!="plain") 
    V_xi_n<-(matrix(out$V_xi_n,ncol=P_n))
  V_xi_e<-out$V_xi_e
  V_z<-t(matrix(out$V_z,ncol=N))
  V_sigma2<-out$V_sigma2
  V_eta<-t(matrix(out$V_eta,ncol=G))
  V_omega2<-out$V_omega2
  V_lambda<-t(matrix(out$V_lambda,ncol=G))
  V_nu<-out$V_nu
  V_alpha<-out$V_alpha
  V_psi2_n<-out$V_psi2_n
  V_psi2_e<-out$V_psi2_e
  
  V_eta<-t(apply(V_eta, 1, "-", apply(V_z, 2, mean)))
  V_z<-t(apply(V_z, 1, "-", apply(V_z, 2, mean)))
  
  variational.params<-list()
  variational.params$net<-variational.start$net
  P_n->variational.params$P_n
  P_e->variational.params$P_e
  model->variational.params$model
  d->variational.params$d
  N->variational.params$N
  NE->variational.params$NE
  NnonE->variational.params$NnonE
  NM->variational.params$NM
  G->variational.params$G
  Y->variational.params$Y
  E->variational.params$E
  nonE->variational.params$nonE
  M->variational.params$M
  numedges->variational.params$numedges
  EnonE->variational.params$EnonE
  diam->variational.params$diam
  hopslist->variational.params$hopslist
  XX_e->variational.params$XX_e
  V_xi_n->variational.params$V_xi_n
  V_xi_e->variational.params$V_xi_e
  V_psi2_n->variational.params$V_psi2_n
  V_psi2_e->variational.params$V_psi2_e
  V_z->variational.params$V_z
  V_sigma2->variational.params$V_sigma2
  V_eta->variational.params$V_eta
  V_lambda->variational.params$V_lambda
  V_omega2->variational.params$V_omega2
  V_nu->variational.params$V_nu
  V_alpha->variational.params$V_alpha
  xi->variational.params$xi
  psi2->variational.params$psi2
  sigma02->variational.params$sigma02
  omega2->variational.params$omega2
  nu->variational.params$nu
  alpha->variational.params$alpha
  inv_sigma02->variational.params$inv_sigma02
  NC->variational.params$NC
  as.logical(out$conv)->variational.params$conv
  seed->variational.params$seed # this is the value the RNG is now using, not the original seed value
  BIC<-vblpcmbic(variational.params)
  BIC->variational.params$BIC
  class(variational.params)<-"vblpcm"
  vblpcmKL(variational.params)
  return(variational.params)
  }
