# reml.Rao.Yu - restricted maximum likelihood estimation of the 
# multivariate Rao-Yu model

reml.Rao.Yu <- function(y, X, M, T, NV=1, vcov_e, maxiter=100, iter.tol=.1e-5, 
         sig2_u = 1, sig2_v=1, rho=.8, rho_u =.4, delta=NULL,
         rho.fixed=NULL, y.include=NULL, ids=NULL, contrast.matrix=NULL,    
         baby.steps=TRUE, dampening=.9, iter.history=FALSE,
         sig2.min.factor=.0001, max.rho_u=.98, max.rho=.98,
         tol=.Machine$double.eps, y.rescale=NULL) {

  fit <- list(model="T: Rao-Yu, REML", convergence=TRUE)

  if(!is.null(y.rescale)) {
    if(y.rescale <= 0) 
      stop("y.rescale must be positive")
    y <- y.rescale*y
    vcov_e <- (y.rescale * y.rescale) * vcov_e
    use.y.rescale <- TRUE
  } else {
    use.y.rescale <- FALSE
    y.rescale <- 1
  }

  if(is.null(delta)) starting.delta <-FALSE else 
    if (length(delta)!=2*NV+1+(NV*(NV-1))/2) {
       starting.delta <- FALSE
       warning("delta of wrong length ignored") 
     } else { 
       starting.delta <- TRUE
       sig2_u <- delta[1:NV]
       sig2_v <- delta[c((NV+1):(2*NV))]
       if(use.y.rescale) {
         sig2_u <- (y.rescale * y.rescale) * sig2_u
         sig2_v <- (y.rescale * y.rescale) * sig2_v
       }
       if(is.null(rho.fixed)) {
         rho <- delta[2*NV+1]
       } else {  
         warning("Value of rho.fixed imposed on specified delta")
         rho <- rho.fixed
       }    
       if(NV>1) rho_u <- delta[(2*NV+2):(2*NV+1+(NV*(NV-1))/2)]
     }
   if(!starting.delta){
      if(length(sig2_u)==1)sig2_u<-rep(sig2_u,NV)
      if(length(sig2_v)==1)sig2_v<-rep(sig2_v,NV)
      if(length(rho_u)==1)rho_u<-rep(rho_u,((NV*(NV-1))/2))
      if(!is.null(rho.fixed))rho<-rho.fixed
   }    
     
   if(!is.null(y.include)) {
     if(length(y.include)!=M)stop("The length of y.include must match M")
     M.save<-M
     y.in <-which(y.include==1)
     M<-length(y.in)
     if(M<M.save) {
       X.save<-X
       y.save<-y
       vcov_e.save<-vcov_e
       index.obs <- rep(0,NV*M*T)
       for (i in 1:M) {
          index.obs[(NV*T*(i-1)+1):(NV*T*i)]<-
                           c((NV*T*(y.in[i]-1)+1):(NV*T*y.in[i]))
       }
       y<-y.save[index.obs]
       X<-X.save[index.obs,]
       vcov_e <-vcov_e.save[index.obs,index.obs]
     }
   }
 
  if(is.null(contrast.matrix)) {
    mi.mat <- diag(NV*T)
    n.contr <- NV*T
    contrast.mse <- contrast.g1 <- contrast.g2 <- contrast.g3 <- 
         contrast.fixed.var <- NULL
  } else {
    if(is.vector(contrast.matrix)) contrast.matrix <- matrix(contrast.matrix,
                                   nrow=length(contrast.matrix),ncol=1)
    if(!is.matrix(contrast.matrix)) 
                         stop("contrast.matrix must be a matrix or a vector")
    if(dim(contrast.matrix)[1]!=NV*T) 
        stop("The number of rows for contrast.matrix must agree with NV*T")
    n.contr <- NV*T + dim(contrast.matrix)[2]
    mi.mat <- matrix(0, nrow=NV*T, ncol=n.contr)
    mi.mat[, 1:(NV*T)] <- diag(NV*T)
    mi.mat[,(NV*T+1):n.contr] <- contrast.matrix
    contrast.mse <- contrast.g1 <- contrast.g2 <- contrast.g3  <- 
         contrast.fixed.var <- matrix(0, nrow=M, ncol=dim(contrast.matrix)[2])
  }
  contrast.est <- contrast.fixed.est <- contrast.wt1 <- contrast.wt2 <- 
         matrix(0, nrow=M, ncol=n.contr)
  
 N<-NV*M*T
 len.rho_u <- (NV*(NV-1))/2
 if(!starting.delta) if (NV>1){ delta<- c(sig2_u=sig2_u[1:NV],
           sig2_v=sig2_v[1:NV],rho=rho,rho_u=rho_u[1:len.rho_u]) 
     } else delta<- c(sig2_u=sig2_u, sig2_v=sig2_v,rho=rho)
 len.delta <- 2*NV+1+len.rho_u
 index.rho <- 2*NV+1
 ix.start <- 1:len.delta
 ix.len <- len.delta
 if(!is.null(rho.fixed)) {
   ix.start<-ix.start[-index.rho]
   ix.len<-ix.len-1
 }
 min.vcov <-min(diag(vcov_e))
 delta.min    <- rep(sig2.min.factor * min.vcov ,2*NV)
 delta.min[index.rho] <- 0
 if (NV>1) delta.min[(2*NV+2):len.delta] <- -max.rho_u
 Gamma_v <- matrix(1,nrow=T,ncol=T)
 d.Gamma_v <- matrix(0,nrow=T,ncol=T)
 A <- orth.c(X)
 y.star <- t(A)%*% y
 if(iter.history) {
     delta.hist<-matrix(0,nrow=len.delta,ncol=maxiter)
     llikelihood.hist<-rep(0,times=maxiter)
     adj.hist<-rep(0,times=maxiter)
     adj.hist[1]<-1
     inf.mat.hist<-array(0,dim=c(len.delta,len.delta,maxiter))
     s.hist<-matrix(0,nrow=len.delta,ncol=maxiter)   
     adj.factor.hist<-rep(0,times=maxiter)
     ix.hist<-matrix(0,nrow=len.delta,ncol=maxiter)   
     warning.hist<-matrix(0,nrow=3,ncol=maxiter)
 }
# initial fit
  rho <- delta[index.rho]
  rho_u <- delta[(index.rho+1):len.delta]
  Gamma_b <- matrix(rep(1:T,times=T),nrow=T, ncol=T)
  Gamma_s <- abs(Gamma_b-t(Gamma_b))
  if(rho > 0) {
     Gamma_u <- rho^(Gamma_s)/(1-rho^2)
     d.Gamma_u <- ((Gamma_s * rho^(Gamma_s) * (1-rho^2)/rho) +
               2*rho* rho^(Gamma_s))/((1-rho^2)^2)
  } else {
     Gamma_u <- diag(T)
     d.Gamma_u <- matrix(0,nrow=T,ncol=T)
     for (i in 2:T) {
         d.Gamma_u[i-1,i]<-1
         d.Gamma_u[i,i-1]<-1 
     }
  }
  u_corr <-diag(1,NV)
  if (NV>1) for (i in 2:NV) for (j in 1:(i-1)){ 
     u_corr[i,j] <- u_corr[j,i] <- rho_u[((i-2)*(i-1))/2+j]
     }
  sqrt_u <- sqrt(delta[1:NV])
  u_base <- ((sqrt_u%*%t(sqrt_u)) * u_corr) %x% Gamma_u
  sqrt_v <- sqrt(delta[(NV+1):(2*NV)])
  v_base <- ((sqrt_v%*%t(sqrt_v)) * u_corr) %x% Gamma_v
  V <- vcov_e + 
       diag(M) %x% (u_base + v_base)
  V.inv <- solve(V,tol=tol)
  V.inv.X <- solve(V,X,tol=tol)
  Xt.V.inv.X <- t(X)%*%V.inv.X
  B_Est   <- solve(Xt.V.inv.X, t(X)%*%solve(V,y,tol=tol),  # (6.2.5)
                   tol=tol)
  res <- y - X %*% B_Est
  V.inv.res <- solve(V,res,tol=tol)
  P <- V.inv - V.inv.X %*% solve(Xt.V.inv.X, t(V.inv.X),   # (p. 101)
         tol=tol)
  C.y.star <- t(A) %*% V %*% A
  llike.c <- -.5 *(determinant(C.y.star)$modulus[1] +  
                 t(y.star) %*% solve(C.y.star,y.star))
  conv.count<-0
  rho.out.count<-0
  for (iter in 1:maxiter) {
     if(iter.history) {
        delta.hist[,iter]<-delta
        llikelihood.hist[iter] <- llike.c + (NV*M*T - dim(X)[2]) *
                                                      log(y.rescale) 
     }

#  compute the current information matrix corresponding to delta
     s <- rep(0,times=len.delta)
     inf.mat <- matrix(0,nrow=len.delta,ncol=len.delta)
     V.j <- array(0,dim=c(T*NV,T*NV,len.delta))
     V.j.m <- array(0,dim=c(N,N,len.delta))
     for (i in 1:NV) { 
       c.m <- matrix(0,nrow=NV,ncol=NV)
       c.m[i,i]<-1
       if(delta[i] > 0) {
          seq.j <- c(1:NV)[-i]
          for (j in seq.j) {
            if (j < i) {
               c.m[i,j] <- rho_u[((i-2)*(i-1))/2+j]*
                             sqrt(delta[j]/delta[i])/2
            } else {
               c.m[i,j] <- rho_u[((j-2)*(j-1))/2+i]*
                             sqrt(delta[j]/delta[i])/2
            }
            c.m[j,i] <- c.m[i,j]
          }
       }  
       V.j[,,i]<- c.m %x% Gamma_u
    }
    
    for (i in 1:NV) { 
       c.m <- matrix(0,nrow=NV,ncol=NV)
       c.m[i,i]<-1
       if(delta[i+NV] > 0) {
          seq.j <- c(1:NV)[-i]
          for (j in seq.j) {
            if (j < i) {
               c.m[i,j] <- rho_u[((i-2)*(i-1))/2+j]*
                             sqrt(delta[j+NV]/delta[i+NV])/2
            } else {
               c.m[i,j] <- rho_u[((j-2)*(j-1))/2+i]*
                             sqrt(delta[j+NV]/delta[i+NV])/2
            }
            c.m[j,i] <- c.m[i,j]
          }
       } 
       V.j[,,i+NV]<- c.m %x% Gamma_v
    }

    V.j[,,index.rho]  <- 
                  ((sqrt_u%*%t(sqrt_u)) * u_corr) %x% d.Gamma_u +
                  ((sqrt_v%*%t(sqrt_v)) * u_corr) %x% d.Gamma_v

    if(NV>1) {
      for (i in 2:NV) for (j in 1:(i-1)){ 
         k<- ((i-2)*(i-1))/2+j
         c.m <- matrix(0,nrow=NV,ncol=NV)
         c.m[i,j]<-c.m[j,i]<-1
         V.j[,,index.rho+k]<- 
                  ((sqrt_u%*%t(sqrt_u)) * c.m) %x% Gamma_u +
                  ((sqrt_v%*%t(sqrt_v)) * c.m) %x% Gamma_v
       }
    }

    for (i in 1:len.delta) {
      for (m in 1:M) {
        r1<-(m-1)*T*NV+1
        r2<- m*T*NV
        for (mp in 1:M) {
          r1p<-(mp-1)*T*NV+1
          r2p<- mp*T*NV
          V.j.m[r1:r2,r1p:r2p,i] <- P[r1:r2,r1p:r2p]%*%as.matrix(V.j[,,i]) 
        }
        s[i] <- s[i] + .5 *
               t(V.inv.res[r1:r2]) %*% as.matrix(V.j[,,i]) %*% 
                 V.inv.res[r1:r2]
      }
      s[i] <- s[i] -.5* sum(diag(V.j.m[,,i])) 
    }
    for (i in 1:len.delta) for (j in (1:i)) {
      inf.mat[i,j] <- .5* sum(sapply(1:N,FUN=function(k) {
                                      t(as.vector(V.j.m[k,,i])) %*% 
                                        as.vector(V.j.m[,k,j])} ))   # (6.2.19)
# #   inf.mat[i,j] <- .5* sum(diag(V.j.m[,,i] %*% V.j.m[,,j]))       # (6.2.19)
      if(i!=j)inf.mat[j,i]<-inf.mat[i,j]
      }

     if(iter.history) {
        inf.mat.hist[,,iter]<-inf.mat
        s.hist[,iter]<-s
     }
     if(starting.delta & maxiter == 1) break
     if(iter > 1 & iter == maxiter) break
     if(iter > 10 | 
         (iter > 1 & starting.delta) ) {
        if((max(abs(delta[1:(2*NV)]-last.delta[1:(2*NV)])) < iter.tol*min.vcov &
            max(abs(delta[index.rho:len.delta]-
                   last.delta[index.rho:len.delta])) < iter.tol)| iter==maxiter) {
               if(conv.count>=1) break else conv.count<-conv.count + 1
        } else conv.count<-0
     }
     if(is.null(rho.fixed)) {
       test.inv <-try(solve(inf.mat,tol=tol))
     } else {
       test.inv <-try(solve(inf.mat[ix.start,ix.start],tol=tol))
     }
     if(class(test.inv)=="try-error") {
        save(y,file="bad_y.Rdata")        # debug
        print("Information matrix has become ill-conditioned")
        print(inf.mat)
        print("Current values of parameters (delta)")
        print(delta)
        print(iter)
        break
     }
     last.delta<-delta
     last.llike<-llike.c
     if(is.null(rho.fixed)) {
       delta.delta.base <- solve(inf.mat,s,tol=tol)
     } else {
       delta.delta.base <- rep(0,len.delta)
       delta.delta.base[ix.start] <- 
              solve(inf.mat[ix.start,ix.start],s[ix.start],tol=tol)
     }  
     if(iter<5 & baby.steps) {
        delta.delta.base <- 2^(iter-5)* delta.delta.base 
     } else if(iter==5 & baby.steps) {
        delta.delta.base <- .75*delta.delta.base 
     } 
     low.change<-TRUE
      
# iterations over adj.iter

     for(adj.iter in c(1:20)) {
        adj.ratio<-rep(1,times=len.delta)
        delta.delta<-delta.delta.base
        ix <-ix.start
        for(ix.index in 1:ix.len) {
           if(index.rho %in% ix ) { 
              if (delta.delta[index.rho] >0) {
                 adj.ratio[index.rho] <- min(1, 
                          (max.rho-last.delta[index.rho])/
                                      delta.delta[index.rho])
              }
           }
           if(NV>1) for(i in c((index.rho+1):len.delta)) {
              if(i %in% ix ) {
                 if (delta.delta[i] > 0) { 
                    adj.ratio[i] <- min(1,(max.rho_u-last.delta[i])/
                                        delta.delta[i])
                 }
              }
           }  
           for (i in ix) {
              if (delta.delta[i] <0) {
                 adj.ratio[i] <- min(1,
                               (delta.min[i]-last.delta[i])/delta.delta[i])
              }
           }
           i.out<-which.min(adj.ratio[ix])
           if(adj.ratio[ix[i.out]]<=tol) {
              ix<-ix[-i.out]
              delta.delta[]<-0
              if(ix.index < ix.len) delta.delta[ix]<- (.5^(adj.iter-1))*
                                solve(inf.mat[ix,ix],s[ix],tol=tol)
           } else if (length(ix)==1 | adj.ratio[ix[i.out]] > .5^(iter/4) | 
                      iter%%2==0) {
              if(any(delta.delta>iter.tol))low.change<-FALSE
              delta.delta<-adj.ratio[ix[i.out]]*delta.delta
              if(iter.history) adj.factor.hist[iter]<-adj.ratio[ix[i.out]]
              break # 
           } else {
              ix<-ix[-i.out]
              delta.delta[]<-0
              if(ix.index < ix.len) delta.delta[ix]<- (.5^(adj.iter-1))*
                                solve(inf.mat[ix,ix],s[ix],tol=tol)
           }
        } # end of loop over ix.index  
        if(length(ix)==0) {
           if(iter.history) {
              adj.factor.hist[iter]<-1
              ix.hist[,iter]<-0
              adj.hist[iter]<-adj.iter
           }
           delta.delta.base<-.5*delta.delta.base
           next # next adj.iter
        }
        if(is.null(rho.fixed)) {
          if(index.rho %in% ix) rho.out.count<-0 else 
                              rho.out.count<-rho.out.count+1
        }
        if(iter.history) {
           ix.hist[,iter]<-0
           ix.hist[ix,iter]<-ix
           adj.hist[iter]<-adj.iter
        }
    
        delta<-last.delta+dampening*delta.delta
        rho <- delta[index.rho]
        rho_u <- delta[(index.rho+1):len.delta]
        Gamma_b <- matrix(rep(1:T,times=T),nrow=T, ncol=T)
        Gamma_s <- abs(Gamma_b-t(Gamma_b))
        if(rho > 0) {
            Gamma_u <- rho^(Gamma_s)/(1-rho^2)
            d.Gamma_u <- ((Gamma_s * rho^(Gamma_s) * (1-rho^2)/rho) +
                   2*rho* rho^(Gamma_s))/((1-rho^2)^2)
        } else {
           Gamma_u <- diag(T)
           d.Gamma_u <- matrix(0,nrow=T,ncol=T)
           for (i in 2:T) {
              d.Gamma_u[i-1,i]<-1
              d.Gamma_u[i,i-1]<-1
           }
        }
         
        u_corr <-diag(1,NV)
        if(NV>1) for (i in 2:NV) for (j in 1:(i-1)){ 
           u_corr[i,j] <- u_corr[j,i] <- rho_u[((i-2)*(i-1))/2+j]
        }
        sqrt_u <- sqrt(delta[1:NV])
        u_base <- ((sqrt_u%*%t(sqrt_u)) * u_corr) %x% Gamma_u
        sqrt_v <- sqrt(delta[(NV+1):(2*NV)])
        v_base <- ((sqrt_v%*%t(sqrt_v)) * u_corr) %x% Gamma_v
        V <- vcov_e + 
              diag(M) %x% (u_base+v_base)
        V.inv <- solve(V,tol=tol)
        V.inv.X <- solve(V,X,tol=tol)
        Xt.V.inv.X <- t(X)%*%V.inv.X
        B_Est   <- solve(Xt.V.inv.X, t(X)%*%solve(V,y,tol=tol),  # (6.2.5)
                      tol=tol)    
        res <- y - X %*% B_Est
        V.inv.res <- solve(V,res,tol=tol)
        P <- V.inv - V.inv.X %*% solve(Xt.V.inv.X, t(V.inv.X),   # (p. 101)
            tol=tol)
        C.y.star <- t(A) %*% V %*% A
        llike.c <- -.5 *(determinant(C.y.star)$modulus[1] +  
                    t(y.star) %*% solve(C.y.star,y.star))
        if(llike.c>last.llike) break
        delta.delta.base<-.5*delta.delta.base
     }
     if(adj.iter==20) {
        if(last.llike-tol>llike.c & ! low.change) {
            if(iter.history)warning.hist[1,iter]<-1
            warning("Iterations halted when likelihood gain ceased")
        }
     }
  }
  if(iter==maxiter & maxiter>1) {
     if(iter.history)warning.hist[2,iter]<-2
     warning("Maxiter=",maxiter," reached, check convergence",
       immediate.=TRUE)
  }
  if(iter>1) {
     if(rho.out.count>1){
        if(iter.history)warning.hist[3,iter]<-3 
        warning("Estimation of rho halted")
     }
  }

  llike <- -.5 *(determinant(V)$modulus[1] +  t(res) %*%  
                       solve(V,res,tol=tol))
  Xt.V.inv.X.inv <- solve(Xt.V.inv.X,tol=tol)
    
  M_term_base <- u_base + v_base
  
  for (m in 1:M){
    r1 <- NV*(m-1)*T+1
    r2 <- NV*m*T
    Xi <- X[r1:r2,]
    for (comp in 1:n.contr) {
      mi<-mi.mat[,comp]
      li<-mi %*% Xi
      M_term.t <- t(mi) %*% M_term_base
      contrast.fixed.est[m,comp] <- li %*% B_Est 
      contrast.est[m,comp] <- contrast.fixed.est[m,comp] +
            M_term.t %*% 
            solve(V[r1:r2,r1:r2],(y[r1:r2]-X[r1:r2,]%*%B_Est),tol=tol)
      contrast.wt1[m,comp] <- t(mi) %*% 
            solve(V[r1:r2,r1:r2],t(M_term.t),tol=tol) / sum(mi*mi)
      contrast.wt2[m,comp] <- contrast.wt1[m,comp] + t(mi) %*% 
            (diag(NV*T) - solve(V[r1:r2,r1:r2],M_term_base,tol=tol)) %*%
            (Xi %*% Xt.V.inv.X.inv %*% t(Xi) %*% 
            solve(V[r1:r2,r1:r2],mi,tol=tol))/sum(mi*mi)
    }
  }
  
# begin MSE calculation
  Gi <- u_base + v_base
  V_bar <- solve(inf.mat, tol=tol)
  d.bi <- matrix(0,nrow=NV*T,ncol=len.delta)
  d.V  <- array(0,dim=c(NV*T,NV*T,len.delta))
  g1 <- g2 <- g3 <- fixed.var <- matrix(0,nrow=M,ncol=n.contr)

  for (i in 1:NV) { 
     c.m <- matrix(0,nrow=NV,ncol=NV)
     c.m[i,i]<-1
     if(delta[i] > 0) {
        seq.j <- c(1:NV)[-i]
        for (j in seq.j) {
          if (j < i) {
             c.m[i,j] <- rho_u[((i-2)*(i-1))/2+j]*
                           sqrt(delta[j]/delta[i])/2
          } else {
             c.m[i,j] <- rho_u[((j-2)*(j-1))/2+i]*
                           sqrt(delta[j]/delta[i])/2
          }
          c.m[j,i] <- c.m[i,j]
        }
     }  
     d.V[,,i]<- c.m %x% Gamma_u
  }

  for (i in 1:NV) { 
     c.m <- matrix(0,nrow=NV,ncol=NV)
     c.m[i,i]<-1
     if(delta[i] > 0) {
        seq.j <- c(1:NV)[-i]
        for (j in seq.j) {
          if (j < i) {
             c.m[i,j] <- rho_u[((i-2)*(i-1))/2+j]*
                           sqrt(delta[j+NV]/delta[i+NV])/2
          } else {
             c.m[i,j] <- rho_u[((j-2)*(j-1))/2+i]*
                           sqrt(delta[j+NV]/delta[i+NV])/2
          }
          c.m[j,i] <- c.m[i,j]
        }
     }  
     d.V[,,i+NV]<- c.m %x% Gamma_v
  }

  d.V[,,index.rho] <- ((sqrt_u%*%t(sqrt_u)) * u_corr) %x% d.Gamma_u +
                      ((sqrt_v%*%t(sqrt_v)) * u_corr) %x% d.Gamma_v
             
  if(NV>1) for (i in 2:NV) for (j in 1:(i-1)){ 
      k<- ((i-2)*(i-1))/2+j
      c.m <- matrix(0,nrow=NV,ncol=NV)
      c.m[i,j]<-c.m[j,i]<-1
      d.V[,,index.rho+k]<- 
                 ((sqrt_u%*%t(sqrt_u)) * c.m) %x% Gamma_u +
                 ((sqrt_v%*%t(sqrt_v)) * c.m) %x% Gamma_v
      }
 
  for (m in 1:M){
    r1 <- NV*(m-1)*T+1
    r2 <- NV*m*T
    Xi <- X[r1:r2,]
    Vi <- V[r1:r2,r1:r2]
    vcov_ei <- vcov_e[r1:r2,r1:r2]
    Vi.inv  <- solve(Vi, tol=tol)
    for (comp in c(1:n.contr)) {
      mi<-mi.mat[,comp]
      li<-mi %*% Xi
      bi <- t(mi)%*%Gi%*%Vi.inv                            # (6.2.10)
      di <- li-bi%*%Xi                                     # (6.2.10)
      for (i in 1:len.delta) {
        d.bi[,i] <- t(t(mi)%*%as.matrix(d.V[,,i])%*%Vi.inv - 
                      t(mi)%*%Gi%*%Vi.inv%*%as.matrix(d.V[,,i])%*%Vi.inv)
      }
      g1[m,comp]  <- t(mi)%*%(Gi-Gi%*%Vi.inv%*%Gi)%*%mi
      g2[m,comp]  <- di%*%Xt.V.inv.X.inv%*%t(di)
      g3[m,comp]  <- sum(diag((t(d.bi)%*%Vi%*%d.bi)%*%V_bar))
      fixed.var[m,comp] <- li%*%Xt.V.inv.X.inv%*%t(li)
    }
   }

  if(!is.null(y.include)) {
    y.out <- which(y.include!=1) 
    MM <- length(y.out)
    if(MM > 0){
      merge.c <- function(y1, y2) {
        y <- rep(0,M.save)         
        y[y.in] <- y1
        y[y.out] <- y2
        return(y)
      }
      merge.m <- function(y1, y2) {
        y <- matrix(0, nrow=M.save, ncol=dim(y1)[2])
        for (i in 1:dim(y1)[2]) {
          y[, i] <- merge.c(as.vector(y1[, i]),
                            as.vector(y2[, i]))
        }
        return(y)
      }
      index.non.obs<-rep(0,NV*MM*T)
      g1.non <- g2.non <- zero.non <- matrix(0, nrow=MM, ncol=n.contr)
      contrast.est.non <- contrast.fixed.non <- 
                matrix(0, nrow=MM, ncol=n.contr)
      for (m in 1:MM) {
        index.non.obs[(NV*T*(m-1)+1):(NV*T*m)]<-
             c((NV*T*(y.out[m]-1)+1):(NV*T*y.out[m]))
        Xi <- as.matrix(X.save[(NV*T*(y.out[m]-1)+1):(NV*T*y.out[m]),])
        for (comp in c(1:n.contr)) {
          mi <- mi.mat[,comp]
          li <- mi %*% Xi
          g1.non[m,comp] <- t(mi)%*%Gi%*%mi
          g2.non[m,comp] <- li%*%Xt.V.inv.X.inv%*%t(li) 
          contrast.est.non[m, comp] <- li %*% B_Est
        } 
      }
      contrast.est <- merge.m(contrast.est, contrast.est.non)
      contrast.fixed.est <- merge.m(contrast.fixed.est, contrast.est.non)
      contrast.wt1 <- merge.m(contrast.wt1, zero.non)
      contrast.wt2 <- merge.m(contrast.wt2, zero.non)
      g1 <- merge.m(g1, g1.non)
      g2 <- merge.m(g2, g2.non)
      g3 <- merge.m(g3, zero.non)
      fixed.var <- merge.m(fixed.var, g2.non)
      M <- M.save
    }
  } 
   
  if(n.contr > NV*T) {
    contrast.g1  <- as.matrix(g1[,(NV*T+1):n.contr])
    contrast.g2  <- as.matrix(g2[,(NV*T+1):n.contr])
    contrast.g3  <- as.matrix(g3[,(NV*T+1):n.contr])
    contrast.mse <- contrast.g1+contrast.g2+2*contrast.g3
    contrast.fixed.var <- as.matrix(fixed.var[,(NV*T+1):n.contr])
  }
   if(NV > 1) {
     eblup <- eblup.g1 <- eblup.g2 <- eblup.g3 <- eblup.wt1 <-
              eblup.wt2 <- est.fixed <- est.fixed.var <-
              matrix(0, nrow=M*T, ncol=NV)
     for (nv in 1:NV) {
       for (t in 1:T) {
         eblup[T*(0:(M-1))+t, nv] <-
                      contrast.est[1:M, T*(nv-1) + t]
         est.fixed[T*(0:(M-1))+t, nv] <-
                      contrast.fixed.est[1:M, T*(nv-1) + t]
         eblup.wt1[T*(0:(M-1))+t, nv] <-
                      contrast.wt1[1:M, T*(nv-1) + t]
         eblup.wt2[T*(0:(M-1))+t, nv] <-
                      contrast.wt2[1:M, T*(nv-1) + t]
         eblup.g1[T*(0:(M-1))+t, nv] <-
                      g1[1:M, T*(nv-1) + t]
         eblup.g2[T*(0:(M-1))+t, nv] <-
                      g2[1:M, T*(nv-1) + t]
         eblup.g3[T*(0:(M-1))+t, nv] <-
                      g3[1:M,T*(nv-1) + t]
         est.fixed.var[T*(0:(M-1))+t, nv] <-
                      fixed.var[1:M, T*(nv-1) + t]
       }
     }
   } else {
     eblup <- eblup.g1 <- eblup.g2 <- eblup.g3 <- eblup.wt1 <-
              eblup.wt2 <- est.fixed <- est.fixed.var <- rep(0, times=M*T)
     for (t in 1:T) {
       eblup[T*(0:(M-1))+t] <- contrast.est[1:M, t]
       est.fixed[T*(0:(M-1))+t] <- contrast.fixed.est[1:M, t]
       eblup.wt1[T*(0:(M-1))+t] <- contrast.wt1[1:M, t]
       eblup.wt2[T*(0:(M-1))+t] <- contrast.wt2[1:M, t]
       eblup.g1[T*(0:(M-1))+t] <- g1[1:M, t]
       eblup.g2[T*(0:(M-1))+t] <- g2[1:M, t]
       eblup.g3[T*(0:(M-1))+t] <- g3[1:M, t]
       est.fixed.var[T*(0:(M-1))+t] <- fixed.var[1:M, t]
     }
   }
   eblup.mse <- eblup.g1 + eblup.g2 + 2*eblup.g3
   if(n.contr > NV*T) {
     contrast.fixed.est <- as.matrix(contrast.fixed.est[, (NV*T+1):n.contr])
     contrast.est <- as.matrix(contrast.est[, (NV*T+1):n.contr])
     contrast.wt1 <- as.matrix(contrast.wt1[, (NV*T+1):n.contr]) 
     contrast.wt2 <- as.matrix(contrast.wt2[, (NV*T+1):n.contr])
   } else {
     contrast.fixed.est <- contrast.est <- contrast.wt1 <- contrast.wt2 <- 
       NULL
   }
  llike <- llike - .5*NV*M*T*log(2*pi)
  llike.c <- llike.c - .5*(NV*M*T - dim(X)[2])*log(2*pi) 

    if(NV>1) parm<-c(rho,delta[1:(2*NV)],rho_u,loglikelihood=llike,
              constrained.ll=llike.c,num.iter=iter) else parm <-
                    c(rho,delta[1:(2*NV)],loglikelihood=llike,
              constrained.ll=llike.c,num.iter=iter)

  B_Est <- c(B_Est)
  names(B_Est) <- attr(X,"dimnames")[[2]]
  if(use.y.rescale) {
    parm[2:(1+2*NV)] <- parm[2:(1+2*NV)]/(y.rescale * y.rescale)
    parm["loglikelihood"] <- parm["loglikelihood"] + NV*M*T*log(y.rescale)
    parm["constrained.ll"] <- parm["constrained.ll"] + (NV*M*T - dim(X)[2]) *
                                                      log(y.rescale) 
    B_Est <- B_Est/y.rescale
    eblup <- eblup/y.rescale
    est.fixed <- est.fixed/y.rescale
    delta[1:(2*NV)] <- delta[1:(2*NV)]/(y.rescale * y.rescale)
    eblup.mse <- eblup.mse/(y.rescale * y.rescale)
    est.fixed.var <- est.fixed.var/(y.rescale * y.rescale)
    if(!is.null(contrast.matrix)) {
      contrast.est <- contrast.est/y.rescale
      contrast.mse <- contrast.mse/(y.rescale * y.rescale)
      contrast.g1 <- contrast.g1/(y.rescale * y.rescale)
      contrast.g2 <- contrast.g2/(y.rescale * y.rescale)
      contrast.g3 <- contrast.g3/(y.rescale * y.rescale)
      contrast.fixed.est <- contrast.fixed.est/y.rescale
      contrast.fixed.var <- contrast.fixed.var/(y.rescale * y.rescale)
    }
    inf.mat[1:(2*NV), ] <- (y.rescale * y.rescale) * inf.mat[1:(2*NV), ]
    inf.mat[, 1:(2*NV)] <- (y.rescale * y.rescale) * inf.mat[, 1:(2*NV)]
    Xt.V.inv.X.inv <- Xt.V.inv.X.inv/(y.rescale * y.rescale)
  }
  if(iter.history) {
    length(llikelihood.hist) <- iter
    delta.hist <- delta.hist[,1:iter]
    length(adj.hist)<-iter
    inf.mat.hist<-inf.mat.hist[,,1:iter]
    s.hist<-s.hist[,1:iter]
    ix.hist<-ix.hist[,1:iter]
    adj.factor.hist<-adj.factor.hist[1:iter]
    warning.hist<-warning.hist[,1:iter]
    if(use.y.rescale) {
      delta.hist[1:(2*NV), ] <- delta.hist[1:(2*NV), ]/(y.rescale * y.rescale)
      inf.mat.hist[1:(2*NV), , ] <- (y.rescale * y.rescale) * 
                                          inf.mat.hist[1:(2*NV), , ]
      inf.mat.hist[, 1:(2*NV), ] <- (y.rescale * y.rescale) * 
                                          inf.mat.hist[, 1:(2*NV), ]

      s.hist[1:(2*NV), ] <- (y.rescale * y.rescale)*s.hist[1:(2*NV), ]
    }
    return(list(eblup=eblup, fit=fit, parm=parm, coef=B_Est,  ids=ids,  
                delta=delta, eblup.mse=eblup.mse, eblup.g1=eblup.g1,   
                eblup.g2=eblup.g2, eblup.g3=eblup.g3, est.fixed=est.fixed, 
                est.fixed.var=est.fixed.var, eblup.wt1=eblup.wt1, 
                eblup.wt2=eblup.wt2, contrast.est=contrast.est, 
                contrast.mse=contrast.mse, contrast.g1=contrast.g1, 
                contrast.g2=contrast.g2, contrast.g3=contrast.g3,
                contrast.fixed.est=contrast.fixed.est,
                contrast.fixed.var=contrast.fixed.var,
                contrast.wt1=contrast.wt1, contrast.wt2=contrast.wt2,
                inf.mat=inf.mat, var.coef=Xt.V.inv.X.inv,
                delta.hist=delta.hist,
                llikelihood.hist=llikelihood.hist,
                adj.hist=adj.hist,
                inf.mat.hist=inf.mat.hist,
                s.hist=s.hist,
                ix.hist=ix.hist,
                adj.factor.hist=adj.factor.hist,
                warning.hist=warning.hist))
    } else {

    return(list(eblup=eblup, fit=fit, parm=parm, coef=B_Est, ids=ids, 
                delta=delta, eblup.mse=eblup.mse, eblup.g1=eblup.g1,   
                eblup.g2=eblup.g2, eblup.g3=eblup.g3, est.fixed=est.fixed, 
                est.fixed.var=est.fixed.var, eblup.wt1=eblup.wt1, 
                eblup.wt2=eblup.wt2, contrast.est=contrast.est, 
                contrast.mse=contrast.mse, contrast.g1=contrast.g1, 
                contrast.g2=contrast.g2, contrast.g3=contrast.g3,
                contrast.fixed.est=contrast.fixed.est,
                contrast.fixed.var=contrast.fixed.var,
                contrast.wt1=contrast.wt1, contrast.wt2=contrast.wt2,
                inf.mat=inf.mat, var.coef=Xt.V.inv.X.inv ))
    }

}
