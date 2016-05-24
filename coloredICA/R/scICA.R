scICA <-
function(Xin,M=dim(Xin)[1],Win=diag(M),tol=0.0001,maxit=20,nmaxit=1,unmixing.estimate="eigenvector",n1,n2,nx01=n1,nx02=n2,h){

 p=dim(Xin)[1]
 if (M>p){
  stop("Number of sources must be less or equal than number 
  of variables")
 }

 if (unmixing.estimate != "eigenvector" && unmixing.estimate != "newton"){
  stop("Methods to estimate the unmixing matrix can be 
  'eigenvector' or 'newton' only")
 }
 
 n = ncol(Xin) 
 Xc =t(scale(t(Xin), center = TRUE, scale = FALSE)) 
 svdcovmat=svd(Xc/sqrt(n))
 K = t(svdcovmat$u %*% diag(1/svdcovmat$d))
 K = K[1:M,] 
 Xc = K %*% Xc

 W1 = Win
 omega1=seq(0,2*pi,len=n1) 
 omega2=seq(0,2*pi,len=n2)
 freq=cbind(rep(omega1,n2),rep(omega2,rep(n1,n2)))
 nx0=nx01*nx02
 x01=seq(0,2*pi,len=nx01)
 x02=seq(0,2*pi,len=nx02)
 freq0=cbind(rep(x02,rep(nx01,nx02)),rep(x01,nx02))
 wlik = -Inf

 rm(list=c("Win"))

 if (unmixing.estimate=="eigenvector"){

  # Here we start we part where the unmixing matrix is estimated
  # through the eigenvector method

  ######################################################
  # Initial Step : Initial spectral density estimation #
  ######################################################

  lim=1
  iter=0
  NInv=0
  index1=rep(1:M,rep(M,M)) 
  index2=rep(1:M,M)
  indx=1:n
  Sc=W1%*%Xc    # Initial sources matrix
  sc=array(0,c(n1,n2,M))
  for (m in 1:M){
   sc[,,m]=t(matrix(Sc[m,],n2,n1))
  }

  vec1=1:n1 -1
  vec2=1:n2 -1
  freq1=(2*pi*vec1)/n1
  freq2=(2*pi*vec2)/n2
  im=complex(real=0,imaginary=1)
  dft=array(0,c(n1,n2,M))
  for (j in vec1){
   for (k in vec2){
    a=cbind(vec1)%*%rbind(rep(freq1[j+1],length(vec2)))
    b=cbind(rep(freq2[k+1],length(vec1)))%*%rbind(vec2)
    for (sou in 1:M){
     dft[j+1,k+1,sou]=sum(sc[,,sou]*exp(-im*(a+b)))
    }
   }
  }

  # Periodogram will be a MxM matrix for each frequency

  period=array(0,c(M,M,n))
  cont=1
  for (j in 1:n2){
   for (i in 1:n1){
    period[,,cont]=1/((2*pi)^2*n)*cbind(dft[i,j,])%*%rbind(Conj(dft[i,j,]))
    cont=cont+1
   }
  }
  l_period=log(period)

  # Estimate of the spectral density

  l_fhat_loc=array(0,c(n,3,M))
  l_fhat_nr=array(0,c(n,3,M))

  # Loc polynomial

  for (m in 1:M){
   l_periodg=Re(l_period[m,m,])
   for (i in 1:n){
    l_fhat_loc[i,,m]=locmulti(freq0[i,],l_periodg,n,freq,h)$ahat
   }
  }
  fhat=exp(l_fhat_loc)

  # N-R

  for (m in 1:M){
   l_periodg=Re(l_period[m,m,])
   for (i in 1:n){
    est=l_fhat_loc[i,,m]
    tol_sde=1
    #iter_nr=0
    while(tol_sde>0.000001){
     est_o=est
     est=est_o-solve(hess(est_o,freq0[i,],l_periodg,n,freq,h))%*%grad(est_o,freq0[i,],l_periodg,n,freq,h)
     tol_sde=sum((est-est_o)^2)
     #iter_nr=iter_nr+1
    }
    l_fhat_nr[i,,m]=est
   }
  }
  fhat=exp(l_fhat_nr)

  g=matrix(0,M,n)
  for (m in 1:M){
   g[m,]=fhat[,1,m]
  }

  Xcc=array(0,c(n1,n2,M))
  for (m in 1:M){
   Xcc[,,m]=matrix(Xc[m,],n1,n2)
  }

  dftx=array(0,c(n1,n2,M))
  for (j in vec1){
   for (k in vec2){
    a=cbind(vec1)%*%rbind(rep(freq1[j+1],length(vec2)))
    b=cbind(rep(freq2[k+1],length(vec1)))%*%rbind(vec2)
    for (sou in 1:M){
     dftx[j+1,k+1,sou]=sum(Xcc[,,sou]*exp(-im*(a+b)))
    }
   }
  }

  X.dft=matrix(0,M,n)
  count=1
  for (j in 1:n2){
   for (i in 1:n1){
    X.dft[,count]=dftx[i,j,]
    count=count+1
   }
  }

  tmp=Re(X.dft[index1,indx] * Conj(X.dft[index2,indx]))*1/((2*pi)^2*n)

  #######################################
  # Starting of the iterative algorithm #
  #######################################

  while( lim>tol & iter<maxit & NInv < nmaxit){
   iter = iter+1
 
   #######################################
   # STEP 1 : Unmixing matrix estimation #
   #######################################
 
   taucount=1
   err=1
   orthoerror=1
   W2=W1
   tau=0.5
   eigenval=rep(0,M)
   while (taucount<60 & err>0.00001 & orthoerror>0.00001){           
    for (j in 1:M){
      Gam=0
      if (j>1){
       for (k in 1:(j-1)){
         nu = matrix(W2[k,],M,1)%*%matrix(W2[k,],1,M)
         Gam = Gam + nu
       }
     }
     tmpmat=t(matrix(tmp%*%matrix(1/g[j,],n,1),M,M))
     tmpV=tmpmat+tau*Gam
     eigenv=eigen(tmpV)
     eigenval[j]=eigenv$values[M]
     W2[j,]=eigenv$vectors[,M]
    }
    orthoerror = sum(sum((W2%*%t(W2) - diag(rep(1,M)))^2))
    err = amari_distance(rerow(W1),rerow(W2))
    taucount = taucount+1
    tau = 2*tau
   }

   # Whittle likelihood of previous step
   wlik2 = -1*sum(eigenval)-1*sum(log(g)) + n * log(abs(det(W2)))

   if (wlik<wlik2){
    Wtmp=W1
    wlik=wlik2
   } else
   print(paste("Spatial Color ICA - Iteration ",iter,": current Whittle likelihood(",wlik2,") is smaller than previous one (",wlik,")."))

   lim=err
   print(paste("Spatial Color ICA - Iteration ",iter,": error is equal to ",lim,sep=""))
   W1=W2
   if ((iter == maxit & NInv < nmaxit)){
    print('Spatial Color ICA: iteration reaches to maximum. Start new iteration.');
    W2=matrix(rnorm(M*M),M,M)
    qrdec=qr(W2)
    W2 = qr.Q(qrdec)
    iter = 0
    NInv = NInv +1
   }

   ########################################
   # STEP 2 : Spectral density estimation #
   ########################################

   Sc = W1 %*% Xc
   sc=array(0,c(n1,n2,M))
   for (m in 1:M){
    sc[,,m]=t(matrix(Sc[m,],n2,n1))
   }

   dft=array(0,c(n1,n2,M))
   for (j in vec1){
    for (k in vec2){
     a=cbind(vec1)%*%rbind(rep(freq1[j+1],length(vec2)))
     b=cbind(rep(freq2[k+1],length(vec1)))%*%rbind(vec2)
     for (sou in 1:M){
      dft[j+1,k+1,sou]=sum(sc[,,sou]*exp(-im*(a+b)))
     }
    }
   }

   # Periodogram will be a MxM matrix for each frequency

   period=array(0,c(M,M,n))

   cont=1
   for (j in 1:n2){
    for (i in 1:n1){
     period[,,cont]=1/((2*pi)^2*n)*cbind(dft[i,j,])%*%rbind(Conj(dft[i,j,]))
     cont=cont+1
    }  
   }

   l_period=log(period)

   # Estimate of the spectral density

   l_fhat_loc=array(0,c(nx0,3,M))
   l_fhat_nr=array(0,c(nx0,3,M))

   # Loc polynomial

   for (m in 1:M){
    l_periodg=Re(l_period[m,m,])
    for (i in 1:nx0){
     l_fhat_loc[i,,m]=locmulti(freq0[i,],l_periodg,n,freq,h)$ahat
    }
   }
   fhat=exp(l_fhat_loc) 

   # N-R

   for (m in 1:M){
    l_periodg=Re(l_period[m,m,])
    for (i in 1:nx0){
     est=l_fhat_loc[i,,m]
     tol_sde=1
     #iter_nr=0
     while(tol_sde>0.000001){
      est_o=est
      est=est_o-solve(hess(est_o,freq0[i,],l_periodg,n,freq,h))%*%grad(est_o,freq0[i,],l_periodg,n,freq,h)
      tol_sde=sum((est-est_o)^2)
      #iter_nr=iter_nr+1
     }
     l_fhat_nr[i,,m]=est
    }
   }
   fhat=exp(l_fhat_nr)

   g=matrix(0,M,n)
   for (m in 1:M){
    g[m,]=fhat[,1,m]
   }

   if (NInv ==nmaxit){
    print("Spatial ColorICA: no convergence")
   }
 
  } # end of the iterative part

  if (wlik>wlik2){ 
   W2=Wtmp
   wlik2=wlik
  }

  W1 = W2

 } # Here we end the "if" related to the eigenvector method

 if (unmixing.estimate=="newton"){

  # Here we start we part where the unmixing matrix is estimated
  # through the eigenvector method

  ######################################################
  # Initial Step : Initial spectral density estimation #
  ######################################################

  lim=1
  iter=0
  NInv=0
  lambda=rep(1,M*(M+1)/2)
  index1=rep(1:M,rep(M,M)) 
  index2=rep(1:M,M)
  tempmx1=matrix(index2,M,M)
  tempmx2=matrix(index1,M,M)
  tempmx1=tempmx1[index1<=index2]
  tempmx2=tempmx2[index1<=index2]
  indx_temp=cbind(diag(rep(1,M*(M+1)/2)),matrix(0,M*(M+1)/2,M*(M-1)/2)) 

  Sc=W1%*%Xc    # Initial sources matrix
  sc=array(0,c(n1,n2,M))
  for (m in 1:M){
   sc[,,m]=t(matrix(Sc[m,],n2,n1))
  }

  vec1=1:n1 -1
  vec2=1:n2 -1
  freq1=(2*pi*vec1)/n1
  freq2=(2*pi*vec2)/n2
  im=complex(real=0,imaginary=1)
  dft=array(0,c(n1,n2,M))
  for (j in vec1){
   for (k in vec2){
    a=cbind(vec1)%*%rbind(rep(freq1[j+1],length(vec2)))
    b=cbind(rep(freq2[k+1],length(vec1)))%*%rbind(vec2)
    for (sou in 1:M){
     dft[j+1,k+1,sou]=sum(sc[,,sou]*exp(-im*(a+b)))
    }
   }
  }

  # Periodogram will be a MxM matrix for each frequency

  period=array(0,c(M,M,n))
  cont=1
  for (j in 1:n2){
   for (i in 1:n1){
    period[,,cont]=1/((2*pi)^2*n)*cbind(dft[i,j,])%*%rbind(Conj(dft[i,j,]))
    cont=cont+1
   }
  }
  l_period=log(period)

  # Estimate of the spectral density

  l_fhat_loc=array(0,c(n,3,M))
  l_fhat_nr=array(0,c(n,3,M))

  # Loc polynomial

  for (m in 1:M){
   l_periodg=Re(l_period[m,m,])
   for (i in 1:n){
    l_fhat_loc[i,,m]=locmulti(freq0[i,],l_periodg,n,freq,h)$ahat
   }
  }
  fhat=exp(l_fhat_loc)

  # N-R

  for (m in 1:M){
   l_periodg=Re(l_period[m,m,])
   for (i in 1:n){
    est=l_fhat_loc[i,,m]
    tol_sde=1
    #iter_nr=0
    while(tol_sde>0.000001){
     est_o=est
     est=est_o-solve(hess(est_o,freq0[i,],l_periodg,n,freq,h))%*%grad(est_o,freq0[i,],l_periodg,n,freq,h)
     tol_sde=sum((est-est_o)^2)
     #iter_nr=iter_nr+1
    }
    l_fhat_nr[i,,m]=est
   }
  }
  fhat=exp(l_fhat_nr)

  g=matrix(0,M,n)
  for (m in 1:M){
   g[m,]=fhat[,1,m]
  }

  Xcc=array(0,c(n1,n2,M))
  for (m in 1:M){
   Xcc[,,m]=matrix(Xc[m,],n1,n2)
  }

  dftx=array(0,c(n1,n2,M))
  for (j in vec1){
   for (k in vec2){
    a=cbind(vec1)%*%rbind(rep(freq1[j+1],length(vec2)))
    b=cbind(rep(freq2[k+1],length(vec1)))%*%rbind(vec2)
    for (sou in 1:M){
     dftx[j+1,k+1,sou]=sum(Xcc[,,sou]*exp(-im*(a+b)))
    }
   }
  }

  X.dft=matrix(0,M,n)
  count=1
  for (j in 1:n2){
   for (i in 1:n1){
    X.dft[,count]=dftx[i,j,]
    count=count+1
   }
  }

  #######################################
  # Starting of the iterative algorithm #
  #######################################

  while( lim>tol & iter<maxit & NInv < nmaxit){
   if (iter==0 && NInv>0){ 
    W2=matrix(rnorm(M*M),M,M)
    qrdec=qr(W2)
    W2 = qr.Q(qrdec)
    lambda=rep(1,M*(M+1)/2)
    W1 = W2
   }
   iter = iter+1
 
   #######################################
   # STEP 1 : Unmixing matrix estimation #
   #######################################
 
   X.wdft = W1%*%X.dft
   W1_inv=solve(W1)
 
   # Score function of Loglikelihood function
 
   Score_l = 2*rowSums(Re(Conj(X.wdft[index1,])*X.dft[index2,])/g[index1,])- n * matrix( W1_inv, M^2, 1)

   # Score function of constraint

   c_e=W1%*%t(W1) - diag(rep(1,M))
   c_e=c_e[index1<=index2] # Lower triangular parts including diag

   J_e1 = W1[tempmx1,index2]*indx_temp[tempmx2,index1]
   J_e2 = W1[tempmx2,index2]*indx_temp[tempmx1,index1]
   J_e = J_e1 + J_e2

   # Hessian Matrix of Loglikelihood function

   temp_l = matrix(0,M^2,M^2)
        
   for  (o in 1:M){
    for  (j in 1:M){
     for  (k in 1:M){
      for  (l in 1:M){
       if ( k==o && o==j && j==l){ 
        temp_l[((o-1)*M+j),((k-1)*M+l)] = lambda[(o*(o+1)/2)]
       }
      }
     }
    }
   } 
                
   temp_h = matrix(0,M^2,M^2)
   tmpm = Re(Conj(X.dft[index1,])*X.dft[index2,])
   for (j in 1:M){
    temp_h[((j-1)*M+1):(j*M),((j-1)*M+1):(j*M)]=2*t(matrix(tmpm%*%cbind(g[j,]^(-1)), M,M))
   }

   rm(tmpm)

   Hessian_l = temp_h - temp_l + n*t(W1_inv[index2,index1]*t(W1_inv[index2,index1]))
   H=rbind(cbind(Hessian_l,-t(J_e)),cbind(-J_e,matrix(0,M*(M+1)/2,M*(M+1)/2)))
   rcondh=rcond(H)

   if (iter == maxit){ 
    print('Spatial ColorICA: iteration reaches to maximum. Start new iteration.')
    iter=0
    NInv=NInv+1
   }
   if  (rcondh< 10^(-15)){
    print('Spatial Color ICA: low condition number. Start new iteration.')
    iter = 0
    NInv = NInv +1
   }
   if  (rcondh>= 10^(-15)){
    V_W = rbind(matrix(t(W1), M^2,1),cbind(lambda)) - solve(H)%*%(rbind( Score_l- (t(J_e)%*%lambda), -cbind(c_e)))
    W2 = t(matrix(V_W[1:M^2],M,M)) 
    lambda = V_W[(M^2+1) : (M*(3*M+1)/2)]
    lim = amari_distance(rerow(W2), rerow(W1))
    W1 = W2
    print(paste("Spatial Color ICA - Iteration ",iter,": error is equal to ",lim,sep=""))
   }
 
   ########################################
   # STEP 2 : Spectral density estimation #
   ########################################
 
   Sc = W1 %*% Xc
   sc=array(0,c(n1,n2,M))
   for (m in 1:M){
    sc[,,m]=matrix(Sc[m,],n1,n2)
   }
   dft=array(0,c(n1,n2,M))
   for (j in vec1){
    for (k in vec2){
     a=cbind(vec1)%*%rbind(rep(freq1[j+1],length(vec2)))
     b=cbind(rep(freq2[k+1],length(vec1)))%*%rbind(vec2)
     for (sou in 1:M){
      dft[j+1,k+1,sou]=sum(sc[,,sou]*exp(-im*(a+b)))
     }
    }
   }

   # Periodogram will be a MxM matrix for each frequency

   period=array(0,c(M,M,n))
   cont=1
   for (j in 1:n2){
    for (i in 1:n1){
     period[,,cont]=1/((2*pi)^2*n)*cbind(dft[i,j,])%*%rbind(Conj(dft[i,j,]))
     cont=cont+1
    } 
   }
   l_period=log(period)

   # Estimate of the spectral density

   l_fhat_loc=array(0,c(nx0,3,M))
   l_fhat_nr=array(0,c(nx0,3,M))

   # Loc polynomial

   for (m in 1:M){
    l_periodg=Re(l_period[m,m,])
    for (i in 1:nx0){
     l_fhat_loc[i,,m]=locmulti(freq0[i,],l_periodg,n,freq,h)$ahat
    }
   }
   fhat=exp(l_fhat_loc)

   # N-R

   for (m in 1:M){
    l_periodg=Re(l_period[m,m,])
    for (i in 1:nx0){
     est=l_fhat_loc[i,,m]
     tol_sde=1
     #iter_nr=0
     while(tol_sde>0.000001){
      est_o=est
      est=est_o-solve(hess(est_o,freq0[i,],l_periodg,n,freq,h))%*%grad(est_o,freq0[i,],l_periodg,n,freq,h)
      tol_sde=sum((est-est_o)^2)
      #iter_nr=iter_nr+1
     }
     l_fhat_nr[i,,m]=est
    }
   }
   fhat=exp(l_fhat_nr)

   g=matrix(0,M,n)
   for (m in 1:M){
    g[m,]=fhat[,1,m]
   }

   if (NInv ==nmaxit){
    print("Spatial Color ICA: no convergence")
   }
 
  }

 } # Here we end the "if" related to the newton method

 wt=W1%*%K

 result = new.env()
 result$W = W1
 result$K = K
 result$A = t(wt) %*% solve(wt %*% t(wt))
 #result$A = solve(t(wt) %*% wt) %*% t(wt)
 result$S = wt %*% Xin
 result$X = Xin
 result$iter = iter
 result$NInv = NInv
 result$den = g 
 as.list(result)

}
