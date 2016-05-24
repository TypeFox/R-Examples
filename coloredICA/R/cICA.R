cICA <-
function(Xin,M=dim(Xin)[1],Win=diag(M),tol=0.0001,maxit=20,nmaxit=1,unmixing.estimate="eigenvector",maxnmodels=100){

 p=dim(Xin)[1]
 if (M>p){
  stop("Number of sources must be less or equal than number 
  of variables")
 }
 
 if (unmixing.estimate != "eigenvector" && unmixing.estimate != "newton"){
  stop("Methods to estimate the unmixing matrix can be 
  'eigenvector' or 'newton' only")
 }

 N = ncol(Xin) 
 Xc =t(scale(t(Xin), center = TRUE, scale = FALSE)) 
 svdcovmat=svd(Xc/sqrt(N))
 K = t(svdcovmat$u %*% diag(1/svdcovmat$d))
 K = K[1:M,]
 Xc = K %*% Xc

 W1 = Win
 wlik = -Inf

 rm(list=c("Win"))

 if (unmixing.estimate=="eigenvector"){

  # Here we start the part where the unmixing matrix is estimated
  # through the eigenvector method

  ######################################################
  # Initial Step : Initial spectral density estimation #
  ######################################################

  freqlength=floor(N/2-1)
  freq = 1:freqlength *2*pi/N
  g = matrix(0,M,freqlength)
  X.dft = t(mvfft(t(Xc)))/sqrt(2*pi*N) # Discrete Fourier Transform
  WXc = W1 %*% Xc
  for( j in 1:M){
   fit = ar.yw(WXc[j,],order.max=maxnmodels)
   if (fit$order==0) g[j,]=fit$var.pred/(2*pi)*rep(1,freqlength)
   else g[j,]=(fit$var.pred/(2*pi))/(abs(1- matrix(fit$ar,1,fit$order)%*%exp(-1i*matrix(1:fit$order,fit$order,1)%*%freq))^2)
  }
  rm(list=c("WXc"))

  lim = 1
  iter=0
  NInv = 0
  index1 = as.double(gl(M,M))
  index2 = as.double(gl(M,1,M^2))
  indx=2:(freqlength+1)
  tmp = Re(X.dft[index1,indx] * Conj(X.dft[index2,indx]))

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
   W2 = W1
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
     tmpmat = t(matrix(tmp%*% matrix(1/g[j,],freqlength,1),M,M))
     tmpV = tmpmat+tau*Gam
     eigenv = eigen(tmpV)
     eigenval[j]=eigenv$values[M]
     W2[j,]=eigenv$vectors[,M]
    }
    orthoerror = sum(sum((W2%*%t(W2) - diag(rep(1,M)))^2))
    err = amari_distance(rerow(W1),rerow(W2))
    taucount = taucount+1
    tau = 2*tau
   } 

   # Whittle likelihood of previous step
   wlik2 = -1*sum(eigenval)-1*sum(log(g)) + N * log(abs(det(W2)))
  
   if (wlik<wlik2){
    Wtmp=W1
    wlik=wlik2
   } else
   print(paste("Color ICA - Iteration ",iter,": current Whittle likelihood(",wlik2,") is smaller than previous one (",wlik,")."))

   lim=err
   print(paste("Color ICA - Iteration ",iter,": error is equal to ",lim,sep=""))
   W1=W2
   if ((iter == maxit & NInv < nmaxit)){
    print('Color ICA: iteration reaches to maximum. Start new iteration.')
    W2=matrix(rnorm(M*M),M,M);
    qrdec=qr(W2)
    W2 = qr.Q(qrdec)
    iter = 0
    NInv = NInv +1
   }

   ########################################
   # STEP 2 : Spectral density estimation #
   ########################################

   WXc = W2 %*% Xc
   for( j in 1:M){
    fit = ar.yw(WXc[j,],order.max=maxnmodels)
    if (fit$order==0) g[j,]=fit$var.pred/(2*pi)*rep(1,freqlength)
    else g[j,] = (fit$var.pred/(2*pi))/(abs(1- matrix(fit$ar,1,fit$order)%*%exp(-1i*matrix(1:fit$order,fit$order,1)%*%freq))^2)
   }

   if (NInv ==nmaxit){
    print("Color ICA: no convergence")
   }
 
  } # end of the iterative part

  if (wlik>wlik2){
   W2=Wtmp
   wlik2=wlik
  }
 
 } # Here we end the "if" related to the eigenvector method

 if (unmixing.estimate=="newton"){

  # Here we start the part where the unmixing matrix is estimated
  # through the newton method

  ######################################################
  # Initial Step : Initial spectral density estimation #
  ######################################################

  freqlength=N
  freq = 1:freqlength *2*pi/N
  g = matrix(0,M,freqlength)
  X.dft = t(mvfft(t(Xc)))/sqrt(2*pi*N) # Discrete Fourier Transform
  WXc = W1 %*% Xc
  for( j in 1:M){
   fit = ar.yw(WXc[j,],order.max=maxnmodels)
   if (fit$order==0) g[j,]=fit$var.pred/(2*pi)*rep(1,freqlength)
   else g[j,]=(fit$var.pred/(2*pi))/(abs(1- matrix(fit$ar,1,fit$order)%*%exp(-1i*matrix(1:fit$order,fit$order,1)%*%freq))^2)
  }

  lim = 1
  iter=0
  NInv = 0
  lambda=rep(1,M*(M+1)/2)
  index1=rep(1:M,rep(M,M)) 
  index2=rep(1:M,M)
  tempmx1=matrix(index2,M,M)
  tempmx2=matrix(index1,M,M)
  tempmx1=tempmx1[index1<=index2]
  tempmx2=tempmx2[index1<=index2]
  indx_temp=cbind(diag(rep(1,M*(M+1)/2)),matrix(0,M*(M+1)/2,M*(M-1)/2))

  #######################################
  # Starting of the iterative algorithm #
  #######################################

  while( lim>tol & iter<maxit & NInv < nmaxit){
   iter = iter+1
   if (iter==0 && NInv>0){ 
    W2=matrix(rnorm(M*M),M,M)
    qrdec=qr(W2)
    W2 = qr.Q(qrdec)
    lambda=rep(1,M*(M+1)/2)
    W1 = W2
   }

   #######################################
   # STEP 1 : Unmixing matrix estimation #
   #######################################

   X.wdft = W1%*%X.dft
   W1_inv=solve(W1)

   # Score function of Loglikelihood function

   Score_l = 2*rowSums(Re(Conj(X.wdft[index1,])*X.dft[index2,])/g[index1,])- N * matrix( W1_inv, M^2, 1)

   Score_l = 2*rowSums(Re(Conj(X.wdft[index1,])*X.dft[index2,])/g[index1,])- N * matrix( W1_inv, M^2, 1)

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

   Hessian_l = temp_h - temp_l + N*t(W1_inv[index2,index1]*t(W1_inv[index2,index1]))
   H=rbind(cbind(Hessian_l,-t(J_e)),cbind(-J_e,matrix(0,M*(M+1)/2,M*(M+1)/2)))
   rcondh=rcond(H)

   if (iter == maxit){ 
    print('Color ICA: iteration reaches to maximum. Start new iteration.')
    iter=0
    NInv=NInv+1
   }
   if  (rcondh< 10^(-15)){
    print('Color ICA: low condition number. Start new iteration.')
    iter = 0
    NInv = NInv +1
   }

   if  (rcondh>= 10^(-15)){
    V_W = rbind(matrix(t(W1), M^2,1),cbind(lambda)) - solve(H)%*%(rbind( Score_l- (t(J_e)%*%lambda), -cbind(c_e)))
    W2 = t(matrix(V_W[1:M^2],M,M)) 
    lambda = V_W[(M^2+1) : (M*(3*M+1)/2)]
    lim = amari_distance(rerow(W2), rerow(W1))
    W1 = W2
    print(paste("Color ICA - Iteration ",iter,": error is equal to ",lim,sep=""))
   }

   ########################################
   # STEP 2 : Spectral density estimation #
   ########################################

   WXc = W2 %*% Xc
   for( j in 1:M){
    fit = ar.yw(WXc[j,],order.max=maxnmodels)
    if (fit$order==0) g[j,]=fit$var.pred/(2*pi)*rep(1,freqlength)
    else g[j,] = (fit$var.pred/(2*pi))/(abs(1- matrix(fit$ar,1,fit$order)%*%exp(-1i*matrix(1:fit$order,fit$order,1)%*%freq))^2)
   }

   if (NInv ==nmaxit){
    print("Color ICA: no convergence")
   }
 
  } # end of the iterative part

 } # Here we end the "if" related to the newton method

 wt=W2%*%K

 result = new.env()
 result$W = W2
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
