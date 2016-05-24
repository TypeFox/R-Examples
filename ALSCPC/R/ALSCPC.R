ALS.CPC <-function(alpha,beta,sigma,epsilon,G,nval,D,S){
    orthonormal<- function(B){
      p <- nrow(B) 
      M <- ncol(B)  
      W <- B  
      if (M > 1){for (i in 2:M){C <- c(crossprod(B[,i],W[,1:(i-1)])) /diag(crossprod(W[,1:(i-1)]))
                                W[,i] <- B[,i] - matrix(W[,1:(i-1)],nrow=p) %*% matrix(C,nrow=i-1)}}
      C <- 1/sqrt(diag(crossprod(W)))
      W <- t(t(W) * C)
      
      return(W)
    }
    #==========================================
    QRdecomposition<-function(A){
      res<-vector("list",length=2)
      Q<-orthonormal(A)
      R<-matrix(nrow=dim(A)[2],ncol=dim(A)[2])
      for(i in 1:dim(A)[1]){
        for( j in 1:dim(A)[2]){
          if(i<j){R[i,j]<-Q[i,]%*%A[,j]}
          if(i>j){R[i,j]<-0 }
          S<-mat.or.vec(dim(A)[1], 1)
          if(i==j){
            if(j==1){R[j,j]<-sqrt(t(A[,i])%*%A[,i])}
            else{
              for(r in 1:(j-1)){S<-S+R[r,j]*Q[r,]};R[j,j]<-sqrt(t(A[,j]-S)%*%(A[,j]-S))}
          }}}
      res[[1]]<-Q
      res[[2]]<-R
      return(res)
    }
    #==================================================
    RetractionQR<-function(D,V){
      result<-QRdecomposition(D+V)[[1]]
      return(result)
    }
    #=============================================
    objectfunctionG<-function(G,nval,D,S){
      #S is list of covariance matrices
      values<-numeric(G)
      for(i in 1:G){values[i]<-nval[i]*log(det(diag(diag(t(D)%*%S[[i]]%*%(D)))))}
      return(sum(values))
    }
    #============================================================
    unconsgradfG<-function(G,nval,D,S){
      p=dim(D)[1]
      #W list of covariance matrices and A list of positive diagonal matrices
      components<-vector("list",length=G)
      Z=mat.or.vec(p, p)
      for(g in 1:G){components[[g]]<-nval[g]*2*t(S[[g]])%*%D%*%solve(diag(diag(t(D)%*%S[[g]]%*%(D))))
                    Z <- Z+components[[g]]}
      return(Z)}
    #=================================================================
    sym<-function(M){(M+t(M))/2}
    #===============================================================
    Projection<-function(Z,D){Z-D%*%sym(t(D)%*%Z)}
    #================================================
    gradfG<-function(G,nval,D,S){Projection(unconsgradfG(G,nval,D,S),D)}
    #================================================================
    frobenius.product<-function(A,B){
      m<-dim(A)[1]
      n<-dim(B)[2]
      C<-matrix(nrow=m,ncol=n)
      for(i in 1:m){
        for(j in 1:n){
          C[i,j]<-A[i,j]*B[i,j]
        }
      }
      return(sum(C))
    }
    #=================== armijo step size ===========================
    Armijo<-function(alpha,beta,sigma,G,nval,D,S){
      m<-0
      repeat{lower<-alpha*sigma*(beta^m)*frobenius.product(gradfG(G,nval,D,S),gradfG(G,nval,D,S))
             Re<-RetractionQR(D,-alpha*beta^m *gradfG(G,nval,D,S))
             upper<-objectfunctionG(G,nval,D,S)-objectfunctionG(G,nval,Re,S)
             if(upper>=lower){armijo<-m;break}
             m<-m+1
      }
      return(alpha*beta^armijo)
    }
    #=====================================================
    Respons<-vector("list")
    Respons[[1]]<-D;t<-NULL
    j<-1
    repeat{
      t[j]<-Armijo(alpha,beta,sigma,G,nval,Respons[[j]],S)
      Respons[[j+1]]<-RetractionQR(Respons[[j]],-t[j]*gradfG(G,nval,Respons[[j]],S))
      if(abs(objectfunctionG(G,nval,Respons[[j]],S)-objectfunctionG(G,nval,Respons[[j+1]],S))<epsilon){fin<-j;break}
      j<-j+1
    }
    return(Respons[[fin]])
  }
