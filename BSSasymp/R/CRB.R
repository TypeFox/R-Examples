CRB <- function(sdf, supp=NULL, A=NULL, eps=1e-04,...)
{
  p <- length(sdf)
  dsdf <- function(x,j){ (sdf[[j]](x+eps)-sdf[[j]](x-eps))/(2*eps)}
  
  if(is.null(supp)) supp <- matrix(c(rep(-8,p),rep(8,p)),ncol=2)
  
  if(is.null(A)) A <- diag(p)  
  
  kap <- NULL
  lambda <- NULL
  for(j in 1:p){ 
    kap[j] <- integrate(Vectorize(function(x){dsdf(x,j)^2/sdf[[j]](x)}),supp[j,1],supp[j,2],...)$value
    lambda[j] <- integrate(Vectorize(function(x){dsdf(x,j)^2*x^2/sdf[[j]](x)}),supp[j,1],supp[j,2],...)$value-1
  }   

  crlb <- matrix(0,p,p)
  for(i in 1:p){
   for(j in 1:p){
    if(i!=j){ 
      crlb[i,j] <- kap[j]/(kap[i]*kap[j]-1)
    }else{
      crlb[i,i] <- 1/(lambda[i])
    }
   } 
  }

  FIM <- matrix(0,p^2,p^2)
  FIM1 <- matrix(0,p,p)
  for(i in 1:p){
   for(j in 1:p){
    if(i==j){
      FIM1 <- tcrossprod(A[,j],A[,i])*lambda[i]
      for(l in 1:p){
       if(l!=i){
        FIM1 <- FIM1+tcrossprod(A[,l],A[,l])*kap[i]
       }
      } 
      FIM[((i-1)*p+1):(i*p),((j-1)*p+1):(j*p)]<-FIM1 
    }else{
     FIM[((i-1)*p+1):(i*p),((j-1)*p+1):(j*p)]<-tcrossprod(A[,j],A[,i])
    }
   }
  }

  EMD <- sum(crlb-diag(crlb)*as.vector(diag(p)))
      
  list(CRLB=crlb, FIM=FIM, EMD=EMD)
}




