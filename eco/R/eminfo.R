##logit and invlogit transformation
logit <- function(X) 
  { 
    Y<-log(X/(1-X))
    Y
   }

invlogit <-function(Y) 
   {
     X<-exp(Y)/(1+exp(Y))
     X
   }

####assuming theta.em
##2 d: mu1, mu2, sig1, sig2, r12
##3 d: mu3, mu1, mu2, sig3, sig1, sig2, r13, r23, r12
param.pack<-function(theta.em, fix.rho=FALSE,r12=0, dim) 
  {
    mu<-rep(0, dim)
    Sig<-matrix(0,dim, dim)

    mu<-theta.em[1:dim]
        
    for (i in 1:dim) 
      Sig[i,i]<-theta.em[dim+i]

   if (!fix.rho) {
      Sig[1,2]<-Sig[2,1]<-theta.em[2*dim+1]*sqrt(Sig[1,1]*Sig[2,2])
      if (dim==3) {
        Sig[1,3]<-Sig[3,1]<-theta.em[2*dim+2]*sqrt(Sig[1,1]*Sig[3,3])
        Sig[2,3]<-Sig[3,2]<-theta.em[2*dim+3]*sqrt(Sig[2,2]*Sig[3,3])
      }
   }
  if (fix.rho) {
    if (dim==2)
       Sig[1,2]<-Sig[2,1]<-r12*sqrt(Sig[1,1]*Sig[2,2])
    if (dim==3) {
      Sig[1,2]<-Sig[2,1]<-theta.em[2*dim+1]*sqrt(Sig[1,1]*Sig[2,2])
      Sig[1,3]<-Sig[3,1]<-theta.em[2*dim+2]*sqrt(Sig[1,1]*Sig[3,3])
      Sig[2,3]<-Sig[3,2]<-r12*sqrt(Sig[2,2]*Sig[3,3])
   }
  }
  return(list(mu=mu, Sigma=Sig))
}

## transformation of BVN parameter into
## Fisher scale or unit scale 
## in 2 D, mu1, mu2, sigma1, sigma2, r12
## in 3 D, mu3, mu1, mu2, sigma3, sigma1, sigma2, sigma31, sigma32, sigma12
param.trans <-function(X, transformation="Fisher") {
  p<-length(X) 
  Y<-rep(0,p)

  if (transformation=="Fisher") {
    if (p<=5) {
      Y[1:2]<-X[1:2]
      Y[3:4]<-log(X[3:4])
      if (p==5) 
        Y[5]<-0.5*log((1+X[5])/(1-X[5]))
     }

     if (p>5) {
       Y[1:3]<-X[1:3]
       Y[4:6]<-log(X[4:6])
       Y[7:8]<-0.5*log((1+X[7:8])/(1-X[7:8]))
       if (p==9)
         Y[9]<-0.5*log((1+X[9])/(1-X[9]))
      }
   }

   if (transformation=="unitscale") {
     if (p<=5) {
    Y[1:2] <- invlogit(X[1:2])
        Y[3:4] <- X[3:4]*exp(2*X[1:2])/(1+exp(X[1:2]))^4
        if (p==5) 
          Y[5] <- X[5]
    }

    if (p>5) {
    Y[1:3]<-invlogit(X[1:3])
        Y[4:6]<-X[4:6]*exp(2*X[4:6])/(1+exp(X[4:6]))^4
        Y[7:8]<-X[7:8]
        if (p==9)
           Y[9]<-X[9]
      }
  }
  return(Y)
}

vec<-function(mat) {
  v<-as.vector(mat, mode="any")
  v
}

tr<-function(mat) {
  trace<-sum(diag(mat))
  trace
}
## I_{com} 
## the gradient function for multivariate normal function
#du.theta and dSig.theta are the first derivative of mu and Sigma 
#with respect to theta
#du.theta[n.u, n.theta]
#dSig.theta[n.u, n.u, n.theta]

d1st.mvn<-function(mu,Sigma, fix.rho=FALSE) {
   #r12, r13,r23 are internal here, 
   # r12 doesn't correspond to cor(w1, w2) in 3d case (intead, r12=>cor(W1,x)
   d<-length(mu)
   p<-d+d+d*(d-1)/2
   u1<-mu[1]
   u2<-mu[2]
   s1<-Sigma[1,1]
   s2<-Sigma[2,2]
   r12<-Sigma[1,2]/sqrt(s1*s2)

   if (d==3) {
   u3<-mu[3]
   s3<-Sigma[3,3]
   r13<-Sigma[1,3]/sqrt(s1*s3)
   r23<-Sigma[2,3]/sqrt(s2*s3)
   }
 
   if (fix.rho) p<-p-1
 
    du.theta<-matrix(0,d,p)
    for (j in 1:d) 
      du.theta[j,j]<-1

    dSig.theta<-array(0,c(d,d,p))
 
      for (i in 1:d) 
        dSig.theta[i,i,d+i]<-1

      dSig.theta[1,2,d+1]<-dSig.theta[2,1,d+1]<-1/2*s1^(-1/2)*s2^(1/2)*r12
      dSig.theta[1,2,d+2]<-dSig.theta[2,1,d+2]<-1/2*s2^(-1/2)*s1^(1/2)*r12
      if (d==3) {
        dSig.theta[1,3,d+1]<-dSig.theta[3,1,d+1]<-1/2*s1^(-1/2)*s3^(1/2)*r13
        dSig.theta[1,3,d+3]<-dSig.theta[3,1,d+3]<-1/2*s3^(-1/2)*s1^(1/2)*r13
        dSig.theta[2,3,d+2]<-dSig.theta[3,2,d+2]<-1/2*s2^(-1/2)*s3^(1/2)*r23
        dSig.theta[2,3,d+3]<-dSig.theta[3,2,d+3]<-1/2*s3^(-1/2)*s2^(1/2)*r23
      }
    
    if (!fix.rho) {
        dSig.theta[1,2,2*d+1]<-dSig.theta[2,1,2*d+1]<-sqrt(s1*s2)
        if (d==3) {
          dSig.theta[1,3,2*d+2]<-dSig.theta[3,1,2*d+2]<-sqrt(s1*s3)
          dSig.theta[2,3,2*d+3]<-dSig.theta[3,2,2*d+3]<-sqrt(s2*s3)
        }
     }
     if (fix.rho) {
        if (d==3) {
          dSig.theta[1,3,2*d+1]<-dSig.theta[3,1,2*d+1]<-sqrt(s1*s3)
          dSig.theta[2,3,2*d+2]<-dSig.theta[3,2,2*d+2]<-sqrt(s2*s3)
        }
     }
    
   return(list(du.theta=du.theta, dSig.theta=dSig.theta))
}

d2nd.mvn<-function(mu,Sigma,  fix.rho=FALSE) {
   #r12, r13,r23 are internal here, 
   # r12 doesn't correspond to cor(w1, w2) in 3d case (intead, r12=>cor(W1,x)
   d<-length(mu)
   p<-d+d+d*(d-1)/2
   u1<-mu[1]
   u2<-mu[2]
   s1<-Sigma[1,1]
   s2<-Sigma[2,2]
   r12<-Sigma[1,2]/sqrt(s1*s2)
   if (d==3) {
   u3<-mu[3]
   s3<-Sigma[3,3]
   r13<-Sigma[1,3]/sqrt(s1*s3)
   r23<-Sigma[2,3]/sqrt(s2*s3)
   }

   if (fix.rho) p<-p-1

   ddu.theta<-array(0,c(d,p,p))

   ddSig.theta<-array(0,c(d,d,p,p))

     ddSig.theta[1,2,d+1,d+1]<-ddSig.theta[2,1,d+1,d+1]<- -1/4*s1^(-3/2)*s2^(1/2)*r12
     ddSig.theta[1,2,d+1,d+2]<-ddSig.theta[2,1,d+1,d+2]<- 1/4*s1^(-1/2)*s2^(-1/2)*r12
     ddSig.theta[1,2,d+2,d+2]<-ddSig.theta[2,1,d+2,d+2]<- -1/4*s1^(1/2)*s2^(-3/2)*r12
 
     if (d==3) {
     ddSig.theta[1,3,d+1,d+1]<-ddSig.theta[3,1,d+1,d+1]<- -1/4*s1^(-3/2)*s3^(1/2)*r13
     ddSig.theta[1,3,d+1,d+3]<-ddSig.theta[3,1,d+1,d+3]<- 1/4*s1^(-1/2)*s3^(-1/2)*r13

     ddSig.theta[2,3,d+2,d+2]<-ddSig.theta[3,2,d+2,d+2]<- -1/4*s2^(-3/2)*s3^(1/2)*r23
     ddSig.theta[2,3,d+2,d+3]<-ddSig.theta[3,2,d+2,d+3]<- 1/4*s2^(-1/2)*s3^(-1/2)*r23

     ddSig.theta[1,3,d+3,d+3]<-ddSig.theta[3,1,d+3,d+3]<- -1/4*s1^(1/2)*s3^(-3/2)*r13
     ddSig.theta[2,3,d+3,d+3]<-ddSig.theta[3,2,d+3,d+3]<- -1/4*s2^(1/2)*s3^(-3/2)*r23

     }

     if (!fix.rho) {
       ddSig.theta[1,2,d+1,2*d+1]<-ddSig.theta[2,1,d+1,2*d+1]<- 1/2*s1^(-1/2)*s2^(1/2)
       ddSig.theta[1,2,d+2,2*d+1]<-ddSig.theta[2,1,d+2,2*d+1]<- 1/2*s1^(1/2)*s2^(-1/2)

       if (d==3) {
         ddSig.theta[1,3,d+1,2*d+2]<-ddSig.theta[3,1,d+1,2*d+2]<- 1/2*s1^(-1/2)*s3^(1/2)
         ddSig.theta[2,3,d+2,2*d+3]<-ddSig.theta[3,2,d+2,2*d+3]<- 1/2*s2^(-1/2)*s3^(1/2)
         ddSig.theta[1,3,d+3,2*d+2]<-ddSig.theta[3,1,d+3,2*d+2]<- 1/2*s1^(1/2)*s3^(-1/2)
         ddSig.theta[2,3,d+3,2*d+3]<-ddSig.theta[3,2,d+3,2*d+3]<- 1/2*s2^(1/2)*s3^(-1/2)
     }
   }
    if (fix.rho) {
       if (d==3) {

       ddSig.theta[1,2,d+1,2*d+1]<-ddSig.theta[2,1,d+1,2*d+1]<- 1/2*s1^(-1/2)*s3^(1/2)
       ddSig.theta[2,3,d+2,2*d+2]<-ddSig.theta[3,2,d+2,2*d+2]<- 1/2*s2^(-1/2)*s3^(1/2)
         ddSig.theta[1,3,d+3,2*d+1]<-ddSig.theta[3,1,d+3,2*d+1]<- 1/2*s1^(1/2)*s3^(-1/2)
         ddSig.theta[2,3,d+3,2*d+2]<-ddSig.theta[3,2,d+3,2*d+2]<- 1/2*s2^(1/2)*s3^(-1/2)
     }
   }

      for (i in 1:(p-1)) 
        for (j in (i+1):p) {
         ddSig.theta[,,j,i]<-ddSig.theta[,,i,j]
         ddu.theta[,j,i]<-ddu.theta[,i,j]
      }

      return(list(ddu.theta=ddu.theta, ddSig.theta=ddSig.theta))
}

##assuming the order of sufficient statistics
## 2d, mean(W1), mean(W2), mean(W1^2) mean(W2^2), mean(W1W2)
## 3d, mean(X), mean(W1), mean(W2), mean(X^2),mean(W1^2) mean(W2^2),
##     mean(XW1), mean(XW2), mean(W1W2)

suff<-function(mu, suff.stat,n) {

   d<-length(mu)
   p<-d+d+d*(d-1)/2 
   u1<-mu[1]
   u2<-mu[2]
   if (d==3)  u3<-mu[3]

   S1<-n*suff.stat[1]
   S2<-n*suff.stat[2]
   S11<-n*suff.stat[d+1]
   S22<-n*suff.stat[d+2]
   S12<-n*suff.stat[2*d+1]
   if (d==3) {
      S3<-n*suff.stat[d]
      S33<-n*suff.stat[2*d]
      S13<-n*suff.stat[2*d+2]
      S23<-n*suff.stat[2*d+3]
   }

   Vv<-rep(0,d)
   Vv[1]<-S1-n*u1
   Vv[2]<-S2-n*u2
   if (d==3) Vv[3]<-S3-n*u3

   Ss<-matrix(0,d,d)
   Ss[1,1]<-S11-2*S1*u1+n*u1^2 
   Ss[2,2]<-S22-2*S2*u2+n*u2^2
   Ss[1,2]<-Ss[2,1]<-S12-S1*u2-S2*u1+n*u1*u2
   if (d==3) {
      Ss[3,3]<-S33-2*S3*u3+n*u3^2
      Ss[1,3]<-Ss[3,1]<-S13-S1*u3-S3*u1+n*u1*u3
      Ss[2,3]<-Ss[3,2]<-S23-S3*u2-S2*u3+n*u2*u3
  }
  return(list(Ss=Ss, Vv=Vv))
}



#du.theta and dSig.theta are the second derivative of mu and Sigma 
#with respect to theta
#ddu.theta[n.u, n.theta, n.theta]
#ddSig.theta[n.u, n.u, n.theta, n.theta]

##comput the gradient vector (expected first derivatives) for MVN
##not actually used here. 

Dcom.mvn<-function(mu, Sigma, suff.stat,n, fix.rho=FALSE) {
  d<-dim(Sigma)[1]
  p<-d*2+0.5*d*(d-1)

  if (fix.rho) { 
    p<-p-1 
  }

  Dcom<-rep(0,p)
  invSigma<-solve(Sigma)

  temp<-suff(mu, suff.stat, n)
  Ss<-temp$Ss
  Vv<-temp$Vv

  temp<-d1st.mvn(mu=mu, Sigma=Sigma, fix.rho=fix.rho)
  du.theta<-temp$du.theta
  dSig.theta<-temp$dSig.theta

  for (i in 1:p)  
   Dcom[i]<- -n/2*t(vec(invSigma))%*%vec(dSig.theta[,,i])+ 0.5*tr(invSigma%*%dSig.theta[,,i]%*%invSigma%*%Ss)+ t(du.theta[,i])%*%invSigma%*%Vv

   Dcom
}


#compute the information matrix of MVN
# -1*second derivatives

Icom.mvn<-function(mu, Sigma, suff.stat,n, fix.rho=FALSE) {
   d<-dim(Sigma)[1]
   p<-d*2+1/2*d*(d-1)

   if (fix.rho) 
     { 
       p<-p-1 
     }

   Icom<-matrix(0,p,p)

   invSigma<-solve(Sigma)

   temp<-suff(mu, suff.stat, n)
   Ss<-temp$Ss
   Vv<-temp$Vv

   temp<-d1st.mvn(mu, Sigma, fix.rho)
   du.theta<-temp$du.theta
   dSig.theta<-temp$dSig.theta

   temp<-d2nd.mvn(mu, Sigma, fix.rho)
   ddu.theta<-temp$ddu.theta
   ddSig.theta<-temp$ddSig.theta

   for (i in 1:p) {
     dinvSig.theta.i<- -invSigma%*%dSig.theta[,,i]%*%invSigma
     for (j in 1:i) {
      dinvSig.theta.j<- -invSigma%*%dSig.theta[,,j]%*%invSigma
      ddinvSig.theta.ij<- -dinvSig.theta.j%*%dSig.theta[,,i]%*%invSigma -invSigma%*%ddSig.theta[,,i,j]%*%invSigma-invSigma%*%dSig.theta[,,i]%*%dinvSig.theta.j
 
       a1<- -n/2*(t(vec(dinvSig.theta.j))%*%vec(dSig.theta[,,i]) + t(vec(invSigma))%*%vec(ddSig.theta[,,i,j]))
                    
       a2<- t(du.theta[,j])%*%dinvSig.theta.i%*%Vv - 0.5*tr(ddinvSig.theta.ij%*%Ss)

     a3<- t(ddu.theta[,i,j])%*%invSigma%*%Vv + t(du.theta[,i])%*%dinvSig.theta.j%*%Vv - n*t(du.theta[,i])%*%invSigma%*%du.theta[,j]

       Icom[i,j]<-a1+a2+a3
    
       if (i!=j) Icom[j,i]<-Icom[i,j]
     }
   }
   -Icom
}

 
###compute the information matrix for various parameter transformation
### "Fisher" transformation (variance stablization?)
### unit scale transformation: first order approximation of mean and var, rho

##express T1 and T2 in more general form

Icom.transform<-function(Icom, Dvec, theta, transformation="Fisher", context, fix.rho)
  {  

      if (!context) {

        mu<-theta[1:2]
        sigma<-theta[3:4]
        rho<-theta[5]
      }
      if (context) {

        mu<-theta[1:3]   # x,w1,w2
        sigma<-theta[4:6] #x, w1, w2
        rho<-theta[7:9]   #r_xw1, r_xw2, r_w1w2
      }
    
    ##T1: d(theta)/d(f(theta)), theta is the MVN parameterization
    ##T2, d2(theta)/d(f(theta))(d(f(theta))')

    ### transformation=Fisher, Icom_normal==>Icom_fisher

    Imat<- -Icom
    n.par<-dim(Imat)[1]

    if (transformation=="Fisher") {
     if (!context) {
         T1<-c(1,1,sigma[1], sigma[2])

         T2<-matrix(0, n.par^2, n.par)
         T2[(2*n.par+3), 3]<-sigma[1]    
         T2[(3*n.par+4), 4]<-sigma[2]     

         if (!fix.rho) {
           T1<-c(T1, (1-(rho[1]^2)))
           T2[(4*n.par+5),5]<- -2*rho[1]*(1-rho[1]^2)
         }
    
        T1<-diag(T1)
     }

     if (context) {
         T1<-c(1,1,1,sigma[1:3],(1-(rho[1:2]^2)))

         T2<-matrix(0, n.par^2, n.par)
         T2[(3*n.par+4), 4]<-sigma[1]    
         T2[(4*n.par+5), 5]<-sigma[2]     
         T2[(5*n.par+6), 6]<-sigma[3]     
         T2[(6*n.par+7),7]<- -2*rho[1]*(1-rho[1]^2)
         T2[(7*n.par+8),8]<- -2*rho[2]*(1-rho[2]^2)

         if (!fix.rho) {
           T1<-c(T1, (1-(rho[3]^2)))
           T2[(8*n.par+9),9]<- -2*rho[3]*(1-rho[3]^2)
         }
    
        T1<-diag(T1)
    }
}
    ### transformation=unitscale, Icom_normal==>Icom_unitscale
   if (transformation=="unitscale") {

      T1<-matrix(0,n.par,n.par)
      T1[1,1]<-exp(-mu[1])*(1+exp(mu[1]))^2
      T1[1,3]<-1/(sigma[1]*2*exp(2*mu[1])*(1+exp(mu[1]))^(-4)*(1-2*(1+exp(mu[1]))^(-1)))
 
      T1[2,2]<-exp(-mu[2])*(1+exp(mu[2]))^2
      T1[2,4]<-1/(sigma[2]*2*exp(2*mu[2])*(1+exp(mu[2]))^(-4)*(1-2*(1+exp(mu[2]))^(-1)))
      
      T1[3,3]<-2*sigma[1]^0.5*(1+exp(mu[1]))^4*exp(-2*mu[1])
      T1[4,4]<-2*sigma[2]^0.5*(1+exp(mu[2]))^4*exp(-2*mu[2])


   #   T2<-matrix(0, n.par^2, n.par)
   #   T2[1,1]<-
   #   T2[(1*n.par+2), (1*n.par+2)]<-


   ##compute T1 and T2 

   }   
 
    Icom.tran<-matrix(NA, n.par, n.par)
    Icom.tran<-T1%*%Imat%*%t(T1)
    
    temp1<-matrix(0,n.par,n.par)
    for (i in 1:n.par)
      for (j in 1:n.par) 
       temp1[i,j]<- Dvec%*%T2[((i-1)*n.par+(1:n.par)),j] 

      Icom.tran<-Icom.tran+temp1     
    return(-Icom.tran)
}


ecoINFO<-function(theta.em, suff.stat, DM, context=TRUE, fix.rho=FALSE, sem=TRUE, r12=0, n)
  {

    if (context) fix.rho<-FALSE
    ndim<-2
    if (context) ndim<-3

    n.var<-2*ndim+ ndim*(ndim-1)/2
 
    n.par<-n.var 
    if (context) {
      n.par<-n.var-2
    }

   if (!context & fix.rho) n.par<-n.par-1

   mu<-param.pack(theta.em, fix.rho=fix.rho, r12=r12,  dim=ndim)$mu
   Sigma<-param.pack(theta.em, fix.rho=fix.rho, r12=r12, dim=ndim)$Sigma
 
  theta.fisher<-param.trans(theta.em)

    Icom<-Icom.mvn(mu=mu, Sigma=Sigma, fix.rho=fix.rho, suff.stat=suff.stat, n=n)
   Dvec<-Dcom.mvn(mu=mu, Sigma=Sigma, fix.rho=fix.rho, suff.stat=suff.stat, n=n)


    theta.icom<-theta.em
    if (fix.rho) theta.icom<-c(theta.em[-n.var], r12)

    Icom.fisher<-Icom.transform(Icom=Icom, Dvec=Dvec, theta=theta.icom, transformation="Fisher", context=context, fix.rho=fix.rho)   


    Vcom.fisher <- solve(Icom.fisher)

      if (!context)  {
      dV <- Vcom.fisher%*%DM%*%solve(diag(1,n.par)-DM)
      Vobs.fisher <- Vcom.fisher+dV }

      ###verify with the parameters.
      ###repartition Icom 
      if (context & !fix.rho) {
       index<-c(1,4,2,3,5,6,7,8,9)
       Itemp<-Icom.fisher[index,index]
       invItemp<-solve(Itemp)
       A1<-invItemp[1:2,1:2]
       A2<-invItemp[1:2,3:9]
       A3<-invItemp[3:9, 1:2]
       A4<-invItemp[3:9, 3:9]
       dV1<-(A4-t(A2)%*%solve(A1)%*%A2)%*%DM%*%solve(diag(rep(1,7))-DM)
       dV<-matrix(0,9,9)
       dV[3:9,3:9]<-dV1
       Vobs.fisher<-invItemp+dV

       index2<-c(1,3,4,2,5,6,7,8,9)
       Vobs.fisher<-Vobs.fisher[index2,index2]
     }

 
 
    Iobs.fisher <- solve(Vobs.fisher)


    ##transform Iobs.fisher to Iobs via delta method
    ##V(theta)=d(fisher^(-1))V(bvn.trans(theta))d(fisher^(-1))'

     if (!context) {
        grad.invfisher <- c(1,1, exp(theta.fisher[3:4]))
        if (! fix.rho)
      grad.invfisher <- c(grad.invfisher,4*exp(2*theta.fisher[5])/(exp(2*theta.fisher[5])+1)^2)
    }

    if (context) {
         grad.invfisher <- c(1,1, 1, exp(theta.fisher[4:6]))
         grad.invfisher <- c(grad.invfisher,4*exp(2*theta.fisher[7:8])/(exp(2*theta.fisher[7:8])+1)^2)
         if (!fix.rho) 
           grad.invfisher <- c(grad.invfisher,4*exp(2*theta.fisher[9])/(exp(2*theta.fisher[9])+1)^2)
    }


    Vobs<-diag(grad.invfisher)%*%Vobs.fisher%*%diag(grad.invfisher)
    Iobs<-solve(Vobs)
    ## obtain a symmetric Cov matrix
    Vobs.sym <- 0.5*(Vobs+t(Vobs))

###unitscale transformation

#theta.unit<-param.trans(theta.em, transformation="unitscale")
#Icom.unit<-Icom.transform(Icom, Dvec,theta.em, transformation="unitscale")
#Vobs.unit<-delta method

if (!context) {
   names(mu)<-c("W1","W2")
   colnames(Sigma)<-rownames(Sigma)<-c("W1","W2")
   names(suff.stat)<-c("S1","S2","S11","S22","S12")
   if (!fix.rho) colnames(DM)<-rownames(DM)<-c("u1","u2","s1","s2","r12")   
   if (fix.rho) colnames(DM)<-rownames(DM)<-c("u1","u2","s1","s2")   
}   
if (context) {
   names(mu)<-c("X","W1","W2")
   colnames(Sigma)<-rownames(Sigma)<-c("X","W1","W2")
   names(suff.stat)<-c("Sx","S1","S2","Sxx","S11","S22","Sx1","Sx2","S12")
   if (!fix.rho) {
    colnames(DM)<-rownames(DM)<-c("u1","u2","s1","s2","r1x","r2x","r12")
    colnames(Icom)<-rownames(Icom)<-c("ux","u1","u2","sx","s1","s2","r1x","r2x","r12")   }  
   if (fix.rho) {
    colnames(DM)<-rownames(DM)<-c("u1","u2","s1","s2","r1x","r2x")   
colnames(Icom)<-rownames(Icom)<-c("ux","u1","u2","sx","s1","s2","r1x","r2x")   }  
}   

colnames(Iobs)<-colnames(Iobs.fisher)<-colnames(Icom.fisher)<-colnames(Vobs)<-colnames(Vobs.sym)<-colnames(Icom)
rownames(Iobs)<-rownames(Iobs.fisher)<-rownames(Icom.fisher)<-rownames(Vobs)<-rownames(Vobs.sym)<-rownames(Icom)

  res.out<-list(mu=mu, Sigma=Sigma, suff.stat=suff.stat, context=context, fix.rho=fix.rho)
  res.out$DM<-DM
    res.out$Icom<-Icom
    res.out$Iobs<-Iobs
    res.out$Fmis<-1-diag(Iobs)/diag(Icom)
    res.out$Vcom<-Vcom<-solve(Icom)
    res.out$Vobs.original<-Vobs
    res.out$VFmis<-1-diag(Vcom)/diag(Vobs)
    res.out$Vobs<-Vobs.sym
    res.out$Icom.trans<-Icom.fisher
    res.out$Iobs.trans<-Iobs.fisher
    res.out$Fmis.trans<-1-diag(Iobs.fisher)/diag(Icom.fisher)
    res.out$Imiss<-res.out$Icom-res.out$Iobs
    res.out$Ieigen<-eigen(res.out$Imiss)[[1]][1]
res.out
}
