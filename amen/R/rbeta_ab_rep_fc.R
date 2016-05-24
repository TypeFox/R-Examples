#' Gibbs sampling of additive row and column effects and regression coefficient
#' with independent replicate relational data
#' 
#' Simulates from the joint full conditional distribution of (a,b,beta),
#' assuming same additive row and column effects and regression coefficient
#' across replicates. 
#' 
#' 
#' @usage rbeta_ab_rep_fc(Z.T,Sab,rho,X.T,s2=1)
#' @param Z.T n x n x T array, with the third dimension for replicates.
#' Each slice of the array is a (latent) normal relational matrix, with
#' multiplicative effects subtracted out
#' @param Sab row and column covariance
#' @param rho dyadic correlation
#' @param X.T n x n x p x T covariate array
#' @param s2 dyadic variance
#' @return \item{beta}{regression coefficients} \item{a}{additive row effects}
#' \item{b}{additive column effects}
#' @author Peter Hoff, Yanjun He
#' @export rbeta_ab_rep_fc
rbeta_ab_rep_fc <-
  function(Z.T,Sab,rho,X.T,s2=1) 
  {
    ###
    N<-dim(X.T)[4]
    p<-dim(X.T)[3] 
    Se<-matrix(c(1,rho,rho,1),2,2)*s2
    iSe2<-mhalf(solve(Se))
    td<-iSe2[1,1] ; to<-iSe2[1,2]
    Sabs<-iSe2%*%Sab%*%iSe2
    tmp<-eigen(Sabs)
    k<-sum(zapsmall(tmp$val)>0 )
    ###
    
    lb<-Qb<-Zr.T<-Zc.T<-Xr.T<-Xc.T<-0
    
    for (t in 1:N){
      Z<-Z.T[,,t]
      X<-array(X.T[,,,t],dim=dim(X.T)[1:3])
      Xr<-apply(X,c(1,3),sum)            # row sum
      Xc<-apply(X,c(2,3),sum)            # col sum
      mX<- apply(X,3,c)                  # design matrix
      mXt<-apply(aperm(X,c(2,1,3)),3,c)  # dyad-transposed design matrix
      XX<-t(mX)%*%mX                     # regression sums of squares
      XXt<-t(mX)%*%mXt
      
      mXs<-td*mX+to*mXt                  # matricized transformed X
      XXs<-(to^2+td^2)*XX + 2*to*td*XXt  # sum of squares for transformed X
      Zs<-td*Z+to*t(Z)
      zr<-rowSums(Zs) ; zc<-colSums(Zs) ; zs<-sum(zc) ; n<-length(zr)
      
      lb<-lb+crossprod(mXs,c(Zs))
      Qb<-Qb+XXs + XX/nrow(mXs)/N
      
      Xsr<-td*Xr + to*Xc  # row sums for transformed X
      Xsc<-td*Xc + to*Xr
      
      Zr.T<-Zr.T+zr
      Zc.T<-Zc.T+zc
      Xr.T<-Xr.T+Xsr
      Xc.T<-Xc.T+Xsc
      
    }
    
    ## dyadic and prior contributions  

    ## row and column reduction
    ab<-matrix(0,nrow(Z),2)
    if(k>0)
    {
      n<-nrow(Z.T[,,1])
      G<-tmp$vec[,1:k] %*% sqrt(diag(tmp$val[1:k],nrow=k))
      K<-matrix(c(0,1,1,0),2,2)
      
      A<-N*n*t(G)%*%G + diag(k)
      B<-N*t(G)%*%K%*%G
      iA0<-solve(A)
      C0<- -solve(A+ n*B)%*%B%*%iA0
      
      iA<-G%*%iA0%*%t(G)
      C<-G%*%C0%*%t(G)
      
      #BigMatrix<-rbind(cbind(n*diag(n),matrix(1,n,n)),cbind(matrix(1,n,n),n*diag(n)))
      #Gn<-G%x%diag(n)
      #V.inv<-N*crossprod(crossprod(BigMatrix,Gn),Gn)+diag(k*n)
      #V<-solve(V.inv)
      #H<-tcrossprod(tcrossprod(Gn,V),Gn)
      H<-iA%x%diag(n)+C%x%matrix(1,nrow=n,ncol=n)
      Hrr<-H[1:n,1:n]
      Hrc<-H[1:n,(n+1):(2*n)]
      Hcr<-H[(n+1):(2*n),1:n]
      Hcc<-H[(n+1):(2*n),(n+1):(2*n)]
      Qb<-Qb-t(Xr.T)%*%Hrr%*%Xr.T-t(Xc.T)%*%Hcr%*%Xr.T-t(Xr.T)%*%Hrc%*%Xc.T-t(Xc.T)%*%Hcc%*%Xc.T
      lb<-lb-t(Xr.T)%*%Hrr%*%Zr.T-t(Xc.T)%*%Hcr%*%Zr.T-t(Xr.T)%*%Hrc%*%Zc.T-t(Xc.T)%*%Hcc%*%Zc.T
      V.b<-solve(Qb)
      M.b<-V.b%*%lb
    }
    ##
    if(dim(X)[3]==0){     beta<-numeric(0) }
    beta<-c(rmvnorm(1,M.b,V.b))
    
    Rr.T<-Zr.T-Xr.T%*%matrix(beta,ncol=1)
    Rc.T<-Zc.T-Xc.T%*%matrix(beta,ncol=1)
    #M.ab<-V%*%(t(G)%x%diag(n))%*%matrix(c(Rr.T,Rc.T),ncol=1)
    #W<-rmvnorm(1,mu=M.ab,V)
    #ab.vec.t<-(G%x%diag(n))%*%matrix(W,ncol=1)
    #ab.vec<-solve(iSe2)%*%matrix(ab.vec.t,nrow=2,byrow=T)
    
    m<-t(t(crossprod(rbind(c(Rr.T),c(Rc.T)),t(iA0%*%t(G)))) + rowSums(sum(Rr.T)*C0%*%t(G)) )
    hiA0<-mhalf(iA0)
    e<-matrix(rnorm(n*k),n,k) 
    w<-m+ t( t(e%*%hiA0) - c(((hiA0-mhalf(iA0+n*C0))/n)%*% colSums(e) ) )
    ab.vec<- t(w%*%t(G)%*%solve(iSe2))
    
    list(beta=beta,a=ab.vec[1,],b=ab.vec[2,] )  
  }


