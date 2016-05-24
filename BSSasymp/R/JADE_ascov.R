ASCOV_JADE <- function(sdf, supp=NULL, A=NULL, ...)
{
  p <- length(sdf)   
  if(is.null(supp)) supp <- matrix(c(rep(-Inf,p),rep(Inf,p)),ncol=2)
  if(is.null(A)) A <- diag(p)
  
  moment3 <- NULL
  moment4 <- NULL
  moment6 <- NULL

  for(j in 1:p){ 
    moment3[j] <- integrate(Vectorize(function(x){sdf[[j]](x)*x^3}),supp[j,1],supp[j,2],...)$value
    moment4[j] <- integrate(Vectorize(function(x){sdf[[j]](x)*x^4}),supp[j,1],supp[j,2],...)$value
    moment6[j] <- integrate(Vectorize(function(x){sdf[[j]](x)*x^6}),supp[j,1],supp[j,2],...)$value   
  }   

  P <- matrix(0,p,p)
  ord <- order(moment4,decreasing=TRUE)
  for(j in 1:p){
    P[j,ord[j]] <- 1
  }  

  moment3 <- moment3[ord]
  moment4 <- moment4[ord]
  moment6 <- moment6[ord]

  kurt <- moment4-3
  
  ASCOV <- matrix(0,p^2,p^2)
  for(i in 1:p){
   for(j in 1:p){
    if(i!=j){ 
      ASVij <- kurt[i]^2*(moment6[i]-moment3[i]^2-kurt[i]^2-6*kurt[i]-9)+kurt[j]^2*(moment6[j]-moment3[j]^2+-6*kurt[j]-9)
 
      ASVij <- ASVij/((kurt[i]^2+kurt[j]^2)^2)
        
      ASCOVij <- -kurt[i]^2*(moment6[i]-moment3[i]^2-kurt[i]^2-6*kurt[i]-9)-kurt[j]^2*(moment6[j]-moment3[j]^2-kurt[j]^2-6*kurt[j]-9-kurt[i]^2)
  
      ASCOVij <- ASCOVij/((kurt[i]^2+kurt[j]^2)^2)
    
      ASCOV <- ASCOV+ASCOVij*kronecker(tcrossprod(diag(p)[,i],diag(p)[,j]),tcrossprod(diag(p)[,j],diag(p)[,i]))+ASVij*kronecker(tcrossprod(diag(p)[,j],diag(p)[,j]),tcrossprod(diag(p)[,i],diag(p)[,i]))

   }  
     
   if(i==j) ASCOV <- ASCOV+0.25*(moment4[j]-1)*kronecker(tcrossprod(diag(p)[,i],diag(p)[,i]),tcrossprod(diag(p)[,i],diag(p)[,i]))    
   } 
  }
  
  EMD <- sum(diag(ASCOV)-diag(ASCOV)*as.vector(diag(p)))  
  W <- crossprod(t(P),solve(A))
  W <- crossprod(diag(sign(rowMeans(W))),W)
  A <- solve(W)
  COV_A <- crossprod(t(tcrossprod(kronecker(diag(p),A),ASCOV)),kronecker(diag(p),t(A)))
  COV_W <- crossprod(t(tcrossprod(kronecker(t(W),diag(p)),ASCOV)),kronecker(W,diag(p)))
  
  list(W=W, COV_W=COV_W, A=A, COV_A=COV_A, EMD=EMD)
}

ASCOV_JADE_est <- function(X, mixed=TRUE)
{
  n <- dim(X)[1]
  p <- dim(X)[2]
  if(mixed){
    W <- JADE(X)$W
  }else W <- diag(p)
  
  X <- tcrossprod(sweep(X,2,colMeans(X)),W)

  moment3 <- NULL 
  moment4 <- NULL
  moment6 <- NULL
 
  for(j in 1:p){
    moment3[j] <- mean(X[,j]^3)
    moment4[j] <- mean(X[,j]^4)
    moment6[j] <- mean(X[,j]^6)
  } 

  P <- matrix(0,p,p)
  ord <- order(moment4,decreasing=TRUE)
  for(j in 1:p){
    P[j,ord[j]] <- 1
  }  
  
  moment3 <- moment3[ord]  
  moment4 <- moment4[ord]
  moment6 <- moment6[ord]
  kurt <- moment4-3

 W <- crossprod(t(P),W)

  ASCOV <- matrix(0,p^2,p^2)
  for(i in 1:p){
   for(j in 1:p){
    if(i!=j){ 
       ASVij <- kurt[i]^2*(moment6[i]-moment3[i]^2-kurt[i]^2-6*kurt[i]-9)+kurt[j]^2*(moment6[j]-moment3[j]^2-6*kurt[j]-9)
 
      ASVij <- ASVij/((kurt[i]^2+kurt[j]^2)^2)
        
      ASCOVij <- -kurt[i]^2*(moment6[i]-moment3[i]^2-kurt[i]^2-6*kurt[i]-9)-kurt[j]^2*(moment6[j]-moment3[j]^2-kurt[j]^2-6*kurt[j]-9-kurt[i]^2)
  
      ASCOVij <- ASCOVij/((kurt[i]^2+kurt[j]^2)^2)
    
      ASCOV <- ASCOV+ASCOVij*kronecker(tcrossprod(diag(p)[,i],diag(p)[,j]),tcrossprod(diag(p)[,j],diag(p)[,i]))+ASVij*kronecker(tcrossprod(diag(p)[,j],diag(p)[,j]),tcrossprod(diag(p)[,i],diag(p)[,i]))

   }  
     
   if(i==j) ASCOV <- ASCOV+0.25*(moment4[j]-1)*kronecker(tcrossprod(diag(p)[,i],diag(p)[,i]),tcrossprod(diag(p)[,i],diag(p)[,i]))         
   } 
  }
  
  A <- solve(W)
  COV_A <- crossprod(t(tcrossprod(kronecker(diag(p),A),ASCOV)),kronecker(diag(p),t(A)))/n  
  COV_W <- crossprod(t(tcrossprod(kronecker(t(W),diag(p)),ASCOV)),kronecker(W,diag(p)))/n
  
  list(W=W, COV_W=COV_W, A=A, COV_A=COV_A)
}

