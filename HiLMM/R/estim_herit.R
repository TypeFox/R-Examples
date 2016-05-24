estim_herit<-
function(Y,W)
{
  n=length(Y)
  N=ncol(W)
  a=n/N
  nb_iter=20
  eta_init=0.5
  eta=rep(0,nb_iter+1)
  eta_init1=0.1
eta_init2=0.5
eta_init3=0.9
eta1=rep(0,nb_iter+1)
eta2=rep(0,nb_iter+1)
eta3=rep(0,nb_iter+1)
eta1[1]=eta_init1
eta2[1]=eta_init2
eta3[1]=eta_init3
  Z=scale(W,center=TRUE,scale=TRUE)
  M=Z%*%t(Z)/N
  O=eigen(M)$vectors
  lambda=eigen(M)$values
  Y_tilde=t(O)%*%Y
   for (nb in 1:nb_iter)
    {A=sum(Y_tilde^2*(lambda-1)/(eta1[nb]*(lambda-1)+1)^2)/sum(Y_tilde^2/(eta1[nb]*(lambda-1)+1))-1/n*sum((lambda-1)/(eta1[nb]*(lambda-1)+1))
     B=((-2*sum(Y_tilde^2*(lambda-1)^2/(eta1[nb]*(lambda-1)+1)^3)*sum(Y_tilde^2/(eta1[nb]*(lambda-1)+1))+(sum(Y_tilde^2*(lambda-1)/(eta1[nb]*(lambda-1)+1)^2))^2)/(sum(Y_tilde^2/(eta1[nb]*(lambda-1)+1)))^2 +1/n*sum((lambda-1)^2/(eta1[nb]*(lambda-1)+1)^2))
     eta1[(nb+1)]=eta1[nb]-A/B
    }
      for (nb in 1:nb_iter)
    {A=sum(Y_tilde^2*(lambda-1)/(eta2[nb]*(lambda-1)+1)^2)/sum(Y_tilde^2/(eta2[nb]*(lambda-1)+1))-1/n*sum((lambda-1)/(eta2[nb]*(lambda-1)+1))
     B=((-2*sum(Y_tilde^2*(lambda-1)^2/(eta2[nb]*(lambda-1)+1)^3)*sum(Y_tilde^2/(eta2[nb]*(lambda-1)+1))+(sum(Y_tilde^2*(lambda-1)/(eta2[nb]*(lambda-1)+1)^2))^2)/(sum(Y_tilde^2/(eta2[nb]*(lambda-1)+1)))^2 +1/n*sum((lambda-1)^2/(eta2[nb]*(lambda-1)+1)^2))
     eta2[(nb+1)]=eta2[nb]-A/B
    }
      for (nb in 1:nb_iter)
    {A=sum(Y_tilde^2*(lambda-1)/(eta3[nb]*(lambda-1)+1)^2)/sum(Y_tilde^2/(eta3[nb]*(lambda-1)+1))-1/n*sum((lambda-1)/(eta3[nb]*(lambda-1)+1))
     B=((-2*sum(Y_tilde^2*(lambda-1)^2/(eta3[nb]*(lambda-1)+1)^3)*sum(Y_tilde^2/(eta3[nb]*(lambda-1)+1))+(sum(Y_tilde^2*(lambda-1)/(eta3[nb]*(lambda-1)+1)^2))^2)/(sum(Y_tilde^2/(eta3[nb]*(lambda-1)+1)))^2 +1/n*sum((lambda-1)^2/(eta3[nb]*(lambda-1)+1)^2))
     eta3[(nb+1)]=eta3[nb]-A/B
    }
      eta_vec=c( eta1[(nb+1)],eta2[(nb+1)],eta3[(nb+1)])
      ind_init=which.min(abs(eta_vec-0.5))
      eta_chap=eta_vec[ind_init]
  w=(lambda-1)/(eta_chap*(lambda-1)+1)
  s_eta=sqrt(2/(sum(w^2)-1/n*sum(w)^2))
  eta_inf=eta_chap -1.96*s_eta
  eta_sup=eta_chap+1.96*s_eta
  list(heritability=eta_chap,CI_low=eta_inf,CI_up=eta_sup,standard_deviation=s_eta)
}
