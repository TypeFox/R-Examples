# lower BIC is preferred
vblpcmbic<-function(v.params)
  {
  N<-v.params$N
  NE<-v.params$NE
  Y<-v.params$Y
  P_n<-v.params$P_n
  P_e<-v.params$P_e
  XX_e<-v.params$XX_e
  V_xi_n<-v.params$V_xi_n
  V_xi_e<-v.params$V_xi_e
  V_psi2_n<-v.params$V_psi2_n
  V_psi2_e<-v.params$V_psi2_e
  V_z<-v.params$V_z
  V_sigma2<-v.params$V_sigma2
  V_omega2<-v.params$V_omega2
  V_eta<-v.params$V_eta
  V_lambda<-v.params$V_lambda
  G<-v.params$G
  d<-v.params$d
  cov1<-XX_e%*%V_xi_e
  cov2<-XX_e%*%V_psi2_e
  if (P_n > 0) 
    {
    tmp<-apply(V_xi_n,1,sum)
    cov1 <- cov1 + c(tmp%*%t(tmp))
    tmp<-rep(sum(V_psi2_n),N) 
    cov2<-cov2+c(tmp%*%t(tmp))
    }
  LL<-0
  for (i in 1:N)
    for (j in (1:N)[-i])
    if (!is.na(Y[i,j]))
      {
      tmp<-cov1[(i-1)*N+j]-sqrt(t(V_z[i,]-V_z[j,])%*%(V_z[i,]-V_z[j,]) + d*(V_sigma2[i]+V_sigma2[j]))
      LL = LL + Y[i,j]*tmp - log(1+exp(tmp + 0.5*cov2[(i-1)*N+j]))
      }
  BIC<-list(
           Y=-2*LL + (N*d+(P_n+P_e-1)*N)*log(N),
	   MBC =
                                if(G>0){
				  tmp=0
				  for (g in 1:G)
				  for (i in 1:N)
				    tmp =
				    tmp-2*sum(V_lambda[g,i]*dnorm(V_z[i,],V_eta[g,],sqrt(V_sigma2[i]+V_omega2[g]),log=TRUE))
				    tmp=tmp+(G+G*d+G-1)*log(N)
				    tmp
                                } else {
                                  -2*sum(dnorm(V_z,0,sqrt(mean(V_z^2)*d),log=TRUE))+1*log(N*d)
                                }
           )
  BIC[["overall"]]<-sum(unlist(BIC))
  BIC
  }
