ASCOV_FastICAdefl <- function(sdf, gs, dgs, Gs=NULL, method="adapt", name=NULL, supp=NULL, A=NULL, ...)
{
  if(method=="adapt"){
   ng <- length(gs)
  }else ng <- 1
  
  if(length(name)!=length(gs)){
   name <- NULL
   for(i in 1:ng){ 
    name[i] <- paste("g",i)
   }
  }
  p <- length(sdf)
  if(is.null(supp)) supp <- matrix(c(rep(-Inf,p),rep(Inf,p)),ncol=2)
  if(is.null(A)) A <- diag(p)
  alpha <- matrix(0,ng,p)
  var_diag <- NULL 

  for(j in 1:p){
    Ex4 <- integrate(Vectorize(function(x){sdf[[j]](x)*x^4}),supp[j,1],supp[j,2])$value
    var_diag[j] <- (Ex4-1)/4

    for(i in 1:ng){
      Eg2 <- integrate(Vectorize(function(x){sdf[[j]](x)*gs[[i]](x)^2}),supp[j,1],supp[j,2],...)$value
      Eg <- integrate(Vectorize(function(x){sdf[[j]](x)*gs[[i]](x)}),supp[j,1],supp[j,2],...)$value
      Egx <- integrate(Vectorize(function(x){sdf[[j]](x)*gs[[i]](x)*x}),supp[j,1],supp[j,2],...)$value
      Edg <- integrate(Vectorize(function(x){sdf[[j]](x)*dgs[[i]](x)}),supp[j,1],supp[j,2],...)$value
     
      alpha[i,j] <- ifelse(abs(Egx-Edg)>1e-06,(Eg2-Eg^2-Egx^2)/(Egx-Edg)^2,Inf)  
    }
  }

  if(method=="adapt"){
   alph <- ifelse(alpha>0,alpha,Inf)
   usedg <- NULL
   for(i in 1:(p-1)){
     mina <- which.min(alph)
     comp <- ceiling(mina/ng)
     gc <- mina-(comp-1)*ng
     usedg[i] <- name[[gc]] 
     alph[,comp] <- Inf
   }  

   ba <- NULL
   for(j in 1:p){
    ba[j] <- min(alpha[,j])
   }
   ord <- order(ba)
   
  }else if(method=="G"){
   Gn <- integrate(Vectorize(function(x){Gs[[1]](x)*dnorm(x)}),-10,10,...)$value  
   EG <- NULL
   for(j in 1:p){
     EG[j] <- integrate(Vectorize(function(x){sdf[[j]](x)*Gs[[1]](x)}),max(-20,supp[j,1]),min(20,supp[j,2]),...)$value
   }   
   ord <- order(abs(EG-Gn),decreasing=TRUE)
   usedg <- rep(name[[1]],p-1)
   ba <- alpha[1,]
  }else{ 
   ord <- 1:p
   usedg <- rep(name[[1]],p-1)
   ba <- alpha[1,]
  }

  P <- diag(p)[ord,]
  bas <- ba[ord]
  var_diags <- var_diag[ord]
 

  ASV <- matrix(0,p,p)
  for(i in 1:p){
   for(j in 1:p){
      if(i<j){ 
        ASV[i,j] <- bas[i]
      }else if(i==j){
        ASV[i,j] <- var_diags[j]
      }else ASV[i,j] <- bas[j]+1 
   }  
  }
  ASCOV <- diag(as.vector(ASV))

  for(i in 1:p){
    for(j in 1:p){
      if(i<j){ 
        ASCOV[(i-1)*p+j,(j-1)*p+i] <- -bas[i]  
       ASCOV[(j-1)*p+i,(i-1)*p+j] <- -bas[i]  
      }else if(i>j){
        ASCOV[(i-1)*p+j,(j-1)*p+i] <- -bas[j]
        ASCOV[(j-1)*p+i,(i-1)*p+j] <- -bas[j]
      }
    }
  }

  EMD <- sum(diag(ASCOV)-diag(ASCOV)*as.vector(diag(p)))
  W <- crossprod(t(P),solve(A))
  W <- crossprod(diag(sign(rowMeans(W))),W)
  A <- solve(W)
  COV_A <- crossprod(t(tcrossprod(kronecker(diag(p),A),ASCOV)),kronecker(diag(p),t(A)))
  COV_W <- crossprod(t(tcrossprod(kronecker(t(W),diag(p)),ASCOV)),kronecker(W,diag(p)))
  
  list(W=W, COV_W=COV_W, A=A, COV_A=COV_A, EMD=EMD, used_gs=usedg)
}

alphas <- function(sdf,gs,dgs,name=NULL,supp=NULL,...)
{ 
  ng<-length(gs)
  if(length(name)!=length(gs)){
   name <- NULL
   for(i in 1:ng){ 
    name[i] <- paste("g",i)
   }
  }
  
  p <- length(sdf)
  cnam<-NULL
  for(j in 1:p){
   cnam[j] <- paste("IC",j)
  }
  
  alpha <- matrix(0,ng,p)
  colnames(alpha) <- cnam
  rownames(alpha) <- name
   
  if(is.null(supp)) supp <- matrix(c(rep(-Inf,p),rep(Inf,p)),ncol=2)

  for(j in 1:p){
    for(i in 1:ng){
      Eg2 <- integrate(Vectorize(function(x){sdf[[j]](x)*gs[[i]](x)^2}),supp[j,1],supp[j,2],...)$value
      Eg <- integrate(Vectorize(function(x){sdf[[j]](x)*gs[[i]](x)}),supp[j,1],supp[j,2],...)$value
      Egx <- integrate(Vectorize(function(x){sdf[[j]](x)*gs[[i]](x)*x}),supp[j,1],supp[j,2],...)$value
      Edg <- integrate(Vectorize(function(x){sdf[[j]](x)*dgs[[i]](x)}),supp[j,1],supp[j,2],...)$value
      
      alpha[i,j] <- ifelse(abs(Egx-Edg)>1e-06,(Eg2-Eg^2-Egx^2)/(Egx-Edg)^2,Inf)
    }
  }
  alpha
}


ASCOV_FastICAdefl_est <- function(X, gs, dgs, Gs=NULL, method="adapt", 
name=NULL, mixed=TRUE)
{
  if(method=="adapt"){
   ng <- length(gs)
  }else ng <- 1

  if(length(name)!=length(gs)){
   name <- NULL
   for(i in 1:ng){ 
    name[i] <- paste("g",i)
   }
  }
  n <- dim(X)[1]
  p <- dim(X)[2] 
  
  if(mixed){
    fI <- adapt_fICA(X,gs,dgs)
    W <- fI$W
    alpha <- fI$alphas
    usedg <- fI$used_gs
  }else{
    alpha <- alphas_data(X,gs,dgs)
    alph <- ifelse(alpha>0,alpha,Inf)
    usedg <- NULL
    for(i in 1:(p-1)){
     mina <- which.min(alph)
     comp <- ceiling(mina/ng)
     gc <- mina-(comp-1)*ng
     usedg[i] <- name[[gc]] 
     alph[,comp] <- Inf
    }  
    W <- diag(p) 
  }  
  
  X <- tcrossprod(sweep(X,2,colMeans(X)),W)   
    
  var_diag <- NULL
  for(j in 1:p){
    var_diag[j] <- (mean(X[,j]^4)-1)/4
  }

  if(method=="adapt"){
   ba <- NULL
   for(j in 1:p){
    ba[j] <- min(alpha[,j])
   }
   ord <- order(ba)
  }else{
   Gn <- integrate(Vectorize(function(x){Gs[[1]](x)*dnorm(x)}),-10,10)$value  
   EG <- NULL
   for(j in 1:p){
     EG[j] <- mean(Gs[[1]](X[,j]))
   }   
   ord <- order(abs(EG-Gn),decreasing=TRUE)
   usedg <- rep(name[[1]],p-1)
   ba <- alpha[1,]
  }

  P <- diag(p)[ord,]    
  bas <- ba[ord]
  var_diags <- var_diag[ord]
   
  W <- P%*%W
 
  ASV <- matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(i<j){ 
        ASV[i,j] <- bas[i]
      }else if(i==j){
        ASV[i,j] <- var_diags[j]
      }else ASV[i,j] <- bas[j]+1 
    }  
  }

  ASCOV <- diag(as.vector(ASV))

  for(i in 1:p){
    for(j in 1:p){
      if(i<j){ 
        ASCOV[(i-1)*p+j,(j-1)*p+i] <- -bas[i]  
        ASCOV[(j-1)*p+i,(i-1)*p+j] <- -bas[i]  
      }else if(i>j){
        ASCOV[(i-1)*p+j,(j-1)*p+i] <- -bas[j]
        ASCOV[(j-1)*p+i,(i-1)*p+j] <- -bas[j]
      }
    }
  }

  A <- solve(W)
  COV_A <- crossprod(t(tcrossprod(kronecker(diag(p),A),ASCOV)),kronecker(diag(p),t(A)))/n  
  COV_W <- crossprod(t(tcrossprod(kronecker(t(W),diag(p)),ASCOV)),kronecker(W,diag(p)))/n
  
  list(W=W, COV_W=COV_W, A=A, COV_A=COV_A, used_gs=usedg)
}


alphas_data <- function(Z, gs, dgs)
{
   ng <- length(gs)
   p <- ncol(Z)
   alpha <- matrix(0,ng,p)
   for(j in 1:p){
    for(i in 1:ng){
     Eg <- mean(gs[[i]](Z[,j]))
     Eg2 <- mean(gs[[i]](Z[,j])^2)
     Egx <- mean(gs[[i]](Z[,j])*Z[,j])
     Edg <- mean(dgs[[i]](Z[,j]))
     alpha[i,j] <- ifelse(abs(Egx-Edg)>1e-06,(Eg2-Eg^2-Egx^2)/(Egx-Edg)^2,Inf)
   } 
   }
   alpha
}


ASCOV_FastICAsym <- function(sdf, G, g, dg, supp=NULL, A=NULL, ...)
{
  p <- length(sdf)
  if(is.null(supp)) supp <- matrix(c(rep(-Inf,p),rep(Inf,p)),ncol=2)
  if(is.null(A)) A <- diag(p)
  var_diag <- NULL
 
  Ex4 <- Eg2 <- Eg <- Egx <- Edg <- EG <- sEGEgx_Edg <- NULL 

  EGn <- integrate(Vectorize(function(x){G(x)*dnorm(x)}),-10,10)$value
  for(j in 1:p){
    Ex4[j] <- integrate(Vectorize(function(x){sdf[[j]](x)*x^4}),supp[j,1],supp[j,2])$value
    var_diag[j] <- (Ex4[j]-1)/4

    Eg2[j] <- integrate(Vectorize(function(x){sdf[[j]](x)*g(x)^2}),supp[j,1],supp[j,2],...)$value
    Eg[j] <- integrate(Vectorize(function(x){sdf[[j]](x)*g(x)}),supp[j,1],supp[j,2],...)$value
    Egx[j] <- integrate(Vectorize(function(x){sdf[[j]](x)*g(x)*x}),supp[j,1],supp[j,2],...)$value
    Edg[j] <- integrate(Vectorize(function(x){sdf[[j]](x)*dg(x)}),supp[j,1],supp[j,2],...)$value
    EG[j] <- integrate(Vectorize(function(x){sdf[[j]](x)*G(x)}),max(-20,supp[j,1]),min(20,supp[j,2]),...)$value-EGn
    sEGEgx_Edg[j] <- sign(EG[j]*(Egx[j]-Edg[j])) 
  } 

  if(sum(sEGEgx_Edg)<sum(abs(sEGEgx_Edg))) stop("at least one bad component")

  ASV <- matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(i!=j){ 
        ASV[i,j] <- (Eg2[i]-Eg[i]^2+Eg2[j]-Eg[j]^2-Egx[i]^2+Edg[j]*(Edg[j]-2*Egx[j]))/((sign(EG[i])*(Egx[i]-Edg[i])+sign(EG[j])*(Egx[j]-Edg[j]))^2)
      }else ASV[i,j] <- var_diag[j]
    }
  }

  ASCOV <- diag(as.vector(ASV))

  for(i in 1:(p-1)){
    for(j in (i+1):p){
       ASCOV[(i-1)*p+j,(j-1)*p+i] <- (-(Eg2[i]-Eg[i]^2)-(Eg2[j]-Eg[j]^2)+Egx[i]^2+
       Egx[j]^2+sign(EG[i])*sign(EG[j])*(Egx[i]-Edg[i])*(Egx[j]-Edg[j]))/
       ((sign(EG[i])*(Egx[i]-Edg[i])+sign(EG[j])*(Egx[j]-Edg[j]))^2)
  
        ASCOV[(j-1)*p+i,(i-1)*p+j] <- ASCOV[(i-1)*p+j,(j-1)*p+i]
    }
  }

  EMD <- sum(diag(ASCOV)-diag(ASCOV)*as.vector(diag(p)))
  W <- solve(A)
  W <- crossprod(diag(sign(rowMeans(W))),W)
  A <- solve(W)
  COV_A <- crossprod(t(tcrossprod(kronecker(diag(p),A),ASCOV)),kronecker(diag(p),t(A)))
  COV_W <- crossprod(t(tcrossprod(kronecker(t(W),diag(p)),ASCOV)),kronecker(W,diag(p)))
  
  list(W=W, COV_W=COV_W, A=A, COV_A=COV_A, EMD=EMD)

}



ASCOV_FastICAsym_est <- function(X, G, g, dg, mixed=TRUE)
{
  n <- dim(X)[1]
  p <- dim(X)[2]
 
  var_diag <- NULL 
  
  if(mixed){
    fI <- fICA(X,g,dg,method="sym")
    W <- fI$W
  }else W <- diag(p)  

  X <- tcrossprod(sweep(X,2,colMeans(X)),W)   
  
  Ex4 <- Eg2 <- Eg <- Egx <- Edg <- EG <- NULL  

  EGn <- integrate(Vectorize(function(x){G(x)*dnorm(x)}),-10,10)$value
  for(j in 1:p){
    Ex4[j] <- mean(X[,j]^4)
    var_diag[j] <- (Ex4[j]-1)/4

    Eg[j] <- mean(g(X[,j]))
    Eg2[j] <- mean(g(X[,j])^2)
    Egx[j] <- mean(g(X[,j])*X[,j])
    Edg[j] <- mean(dg(X[,j]))
    EG[j] <- mean(G(X[,j]))-EGn   
  }

  ASV <- matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(i!=j){ 
        ASV[i,j] <- (Eg2[i]-Eg[i]^2+Eg2[j]-Eg[j]^2-Egx[i]^2+Edg[j]*(Edg[j]-2*Egx[j]))/((sign(EG[i])*(Egx[i]-Edg[i])+sign(EG[j])*(Egx[j]-Edg[j]))^2)
      }else ASV[i,j] <- var_diag[j]
    }
  }

  ASCOV <- diag(as.vector(ASV))

  for(i in 1:(p-1)){
    for(j in (i+1):p){
       ASCOV[(i-1)*p+j,(j-1)*p+i] <- (-(Eg2[i]-Eg[i]^2)-(Eg2[j]-Eg[j]^2)+Egx[i]^2+
       Egx[j]^2+sign(EG[i])*sign(EG[j])*(Egx[i]-Edg[i])*(Egx[j]-Edg[j]))/
       ((sign(EG[i])*(Egx[i]-Edg[i])+sign(EG[j])*(Egx[j]-Edg[j]))^2)
  
        ASCOV[(j-1)*p+i,(i-1)*p+j] <- ASCOV[(i-1)*p+j,(j-1)*p+i]
    }
  }

  A <- solve(W)
  COV_A <- crossprod(t(tcrossprod(kronecker(diag(p),A),ASCOV)),kronecker(diag(p),t(A)))/n  
  COV_W <- crossprod(t(tcrossprod(kronecker(t(W),diag(p)),ASCOV)),kronecker(W,diag(p)))/n
  
  list(W=W, COV_W=COV_W, A=A, COV_A=COV_A)
}


ASCOV_FastICAsym2 <- function(sdf, G, g, dg, supp=NULL, A=NULL, ...)
{
  p <- length(sdf)
  if(is.null(supp)) supp <- matrix(c(rep(-Inf,p),rep(Inf,p)), ncol=2)
  if(is.null(A)) A <- diag(p)
  var_diag <- NULL
 
  Ex4 <- Eg2 <- Eg <- Egx <- Edg <- EG <- sEGEgx_Edg <- NULL  

  EGn <- integrate(Vectorize(function(x){G(x)*dnorm(x)}),-10,10)$value
  for(j in 1:p){
    Ex4[j] <- integrate(Vectorize(function(x){sdf[[j]](x)*x^4}),supp[j,1],supp[j,2])$value
    var_diag[j] <- (Ex4[j]-1)/4

    Eg2[j] <- integrate(Vectorize(function(x){sdf[[j]](x)*g(x)^2}),supp[j,1],supp[j,2],...)$value
    Eg[j] <- integrate(Vectorize(function(x){sdf[[j]](x)*g(x)}),supp[j,1],supp[j,2],...)$value
    Egx[j] <- integrate(Vectorize(function(x){sdf[[j]](x)*g(x)*x}),supp[j,1],supp[j,2],...)$value
    Edg[j] <- integrate(Vectorize(function(x){sdf[[j]](x)*dg(x)}),supp[j,1],supp[j,2],...)$value
    EG[j] <- integrate(Vectorize(function(x){sdf[[j]](x)*G(x)}),max(-20,supp[j,1]),min(20,supp[j,2]),...)$value-EGn
    sEGEgx_Edg[j] <- sign(EG[j]*(Egx[j]-Edg[j])) 
  } 

  if(sum(sEGEgx_Edg)<sum(abs(sEGEgx_Edg))) stop("at least one bad component")

  ASV <- matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(i!=j){ 
        ASV[i,j] <- (EG[i]^2*(Eg2[i]-Eg[i]^2-Egx[i]^2)+EG[j]^2*(Eg2[j]-Eg[j]^2-
2*Egx[j]*Edg[j]+Edg[j]^2))/((EG[i]*(Egx[i]-Edg[i])+EG[j]*(Egx[j]-Edg[j]))^2)
      }else ASV[i,j] <- var_diag[j]
    }
  }

  ASCOV <- diag(as.vector(ASV))

  for(i in 1:(p-1)){
    for(j in (i+1):p){
       ASCOV[(i-1)*p+j,(j-1)*p+i] <- (-EG[i]^2*(Eg2[i]-Eg[i]^2-Egx[i]^2)-EG[j]^2*(Eg2[j]-
Eg[j]^2-Egx[j]^2)+EG[i]*EG[j]*(Egx[i]-Edg[i])*(Egx[j]-Edg[j]))/((EG[i]*(Egx[i]-Edg[i])+EG[j]*(Egx[j]-Edg[j]))^2)
       
        ASCOV[(j-1)*p+i,(i-1)*p+j] <- ASCOV[(i-1)*p+j,(j-1)*p+i]
    }
  }

  EMD <- sum(diag(ASCOV)-diag(ASCOV)*as.vector(diag(p)))
  W <- solve(A)
  W <- crossprod(diag(sign(rowMeans(W))),W)
  A <- solve(W)
  COV_A <- crossprod(t(tcrossprod(kronecker(diag(p),A),ASCOV)),kronecker(diag(p),t(A)))
  COV_W <- crossprod(t(tcrossprod(kronecker(t(W),diag(p)),ASCOV)),kronecker(W,diag(p)))
  
  list(W=W, COV_W=COV_W, A=A, COV_A=COV_A, EMD=EMD)

}



ASCOV_FastICAsym2_est <- function(X, G, g, dg, mixed=TRUE)
{
  n <- dim(X)[1]
  p <- dim(X)[2]
 
  var_diag <- NULL
  
  if(mixed){
    fI <- fICA(X,g,dg,G,method="sym2")
    W <- fI$W
  }else W <- diag(p)  

  X <- tcrossprod(sweep(X,2,colMeans(X)),W)   
  
  Ex4 <- Eg2 <- Eg <- Egx <- Edg <- EG <- NULL

  EGn <- integrate(Vectorize(function(x){G(x)*dnorm(x)}),-10,10)$value
  for(j in 1:p){
    Ex4[j] <- mean(X[,j]^4)
    var_diag[j] <- (Ex4[j]-1)/4

    Eg[j] <- mean(g(X[,j]))
    Eg2[j] <- mean(g(X[,j])^2)
    Egx[j] <- mean(g(X[,j])*X[,j])
    Edg[j] <- mean(dg(X[,j]))
    EG[j] <- mean(G(X[,j]))-EGn   
  }

  ASV <- matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(i!=j){ 
        ASV[i,j] <- (EG[i]^2*(Eg2[i]-Eg[i]^2-Egx[i]^2)+EG[j]^2*(Eg2[j]-Eg[j]^2-2*Egx[j]*Edg[j]+Edg[j]^2))/((EG[i]*(Egx[i]-Edg[i])+EG[j]*(Egx[j]-Edg[j]))^2)
      }else ASV[i,j] <- var_diag[j]
    }
  }

  ASCOV <- diag(as.vector(ASV))

  for(i in 1:(p-1)){
    for(j in (i+1):p){
       ASCOV[(i-1)*p+j,(j-1)*p+i] <- (-EG[i]^2*(Eg2[i]-Eg[i]^2-Egx[i]^2)-EG[j]^2*(Eg2[j]-
Eg[j]^2-Egx[j]^2)+EG[i]*EG[j]*(Egx[i]-Edg[i])*(Egx[j]-Edg[j]))/
((EG[i]*(Egx[i]-Edg[i])+EG[j]*(Egx[j]-Edg[j]))^2)
       
        ASCOV[(j-1)*p+i,(i-1)*p+j] <- ASCOV[(i-1)*p+j,(j-1)*p+i]
    }
  }

  A <- solve(W)
  COV_A <- crossprod(t(tcrossprod(kronecker(diag(p),A),ASCOV)),kronecker(diag(p),t(A)))/n  
  COV_W <- crossprod(t(tcrossprod(kronecker(t(W),diag(p)),ASCOV)),kronecker(W,diag(p)))/n
  
  list(W=W, COV_W=COV_W, A=A, COV_A=COV_A)
}


