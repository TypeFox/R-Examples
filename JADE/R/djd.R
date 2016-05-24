# Function for deflation based JD
#
# input: 
#   X = p,p,k array with k pxp matrices which should be diagonalized
#   G = for optimization, options are "pow" and "log", default is "pow"
#   r = power used, if G="pow", default is 2 
#   eps = convergence criterion
#   maxiter = maximum number of iterations
#
# output:
#   orthogonal matrix

djd <- function(X, G="max", r=2, eps = 1e-06, maxiter = 500)
    {
    G <- match.arg(G, c("max", "pow", "log"))
    
    WN <- eigen(X[,,1])$vectors   
    p <- dim(X)[1]
    k <- dim(X)[3]
    W <- matrix(0,p,p)
    
    W <- switch(G,
        "pow"={
               djd.pow(X=X, W=W, r=r, WN=WN, k=k, p=p, eps=eps, maxiter=maxiter)
               }
        ,
        "log"={
              djd.log(X=X, W=W, WN=WN, k=k, p=p, eps=eps, maxiter=maxiter)
              }
        ,
        "max"={
               djd.max(X=X, W=W, r=r, k=k, p=p, eps=eps, maxiter=maxiter)
               }
        )
    
     W
     } 

########################################################
#
# djd subfunctions
#
djd.pow <- function(X,W,r,WN,k,p,eps,maxiter)
    {
    for (i in 1:p){
                    wn <- (WN[,i, drop=FALSE])
                    for (it in 1:maxiter){
                        w <- wn
                        wn <- matrix(0,p,1)
                        for (mi in 1:k){
                            wXmiw <- as.numeric(crossprod(w, X[,,mi]) %*% w)
                            wn <- wn + 2*sign(wXmiw) * r*(abs(wXmiw))^(r-1) * X[,,mi] %*% w}
                        if(floor(it/5) == it/5){
                         wn <- 0.5*wn + 0.5*w
                        }
                        wn <- wn- crossprod(W)%*%wn
                        wn <- wn / sqrt(sum(wn^2)) 
                        if (sqrt(sum((w-wn)^2))<eps || sqrt(sum((w+wn)^2))<eps) break
                        if (it==maxiter) stop("no convergence reached")
                        }
                    W[i,] <- t(wn)
                    }
    t(W)
    }
    
        
djd.log <- function(X,W,WN,k,p,eps,maxiter)
    {
    for (i in 1:p){
                    wn <- (WN[,i, drop=FALSE])
                    for (it in 1:maxiter){
                        w <- wn
                        wn <- matrix(0,p,1)
                        for (mi in 1:k){
                            wXmi <- crossprod(w, X[,,mi])
                            wn <- wn +  1/as.numeric(wXmi %*% w)* t(wXmi)}
                        wn <- wn- crossprod(W)%*%wn
                        wn <- wn / sqrt(sum(wn^2)) 
                        if (sqrt(sum((w-wn)^2))<eps || sqrt(sum((w+wn)^2))<eps) break
                        if (it==maxiter) stop("no convergence reached")
                        }
                    W[i,] <- t(wn)
                    }
    t(W)
    }


rand_orth <- function(p, W0 = NULL)
{
   W <- matrix(0,p,p)
   if(is.matrix(W0)){
     k <- ncol(W0)
     W[,1:k] <- W0
   }else k <- 0
 
   V <- matrix(rnorm(p*(p-k)),p,p-k)
 
   for(i in 1:(p-k)){
     V[,i] <- V[,i]-tcrossprod(W) %*% V[,i]  
     V[,i] <- V[,i]/sqrt(sum(V[,i]^2))
     W[,k+i] <- V[,i]
   }
 
   W
}


#   nr = number of candidates for the initial value / p
djd.max <- function(X,W,r,k,p,nr=100,eps,maxiter)
    {
    for (i in 1:(p-1)){
        if(i>1){
          U <- rand_orth(p,as.matrix(W[,1:(i-1)]))[,i:p]
        }else U <- diag(p) 
  
        inits <- NULL
        for(j in 1:nr){
          inits <- rbind(inits,rand_orth(p-i+1)) 
        }
        

        winit <- matrix(0,nrow(inits),p)
        for(j in 1:nrow(inits)){
          winit[j,] <- U %*% inits[j,] 
        }
  
        dsum <- NULL
        for(j in 1:nrow(inits)){
          dsum[j] <- 0
          for(l in 1:k){
            dsum[j] <- dsum[j]+abs((t(winit[j,]) %*% X[,,l] %*% winit[j,])^r)
          }
        } 
        wn <- winit[which(dsum==max(dsum)),]
    

                    for (it in 1:maxiter){
                        w <- wn
                        wn <- matrix(0,p,1)
                        for (mi in 1:k){
                            wXmiw <- as.numeric(crossprod(w, X[,,mi]) %*% w)
                            wn <- wn + 2*sign(wXmiw) * r*abs(wXmiw)^(r-1) * X[,,mi] %*% w}
                        if(floor(it/5)==it/5){
                         wn<-0.5*wn+0.5*w
                        }
                        wn <- wn- tcrossprod(W) %*% wn
                        wn <- wn / sqrt(sum(wn^2)) 
                        if (sqrt(sum((w-wn)^2))<eps || sqrt(sum((w+wn)^2))<eps) break
                        if (it==maxiter) stop("no convergence reached")
                        }
                    W[,i] <- t(wn)
                    }
   
         W[,p] <- rand_orth(p,W[,1:(p-1)])[,p]
        
    W
    }


