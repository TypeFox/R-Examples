

mpinv <- function(X)                                                                                                                          
{                                                                                                                                             
    Eps <- 100 * .Machine$double.eps                                                                                                       
                                                                                                                                              
    # singular value decomposition                                                                                                            
                                                                                                                                              
    s <- svd(X)                                                                                                                              
    d <- s$d                                                                                                                                 
    m <- length(d)
                                                                                                                                              
    if (!(is.vector(d)))                                                                                                                      
        return(t(s$v%*%(1/d)%*%t(s$u)))                                                                                                       
                                                                                                                                              
    d <- d[d > Eps]                                                                                                                          
    notnull <- length(d)                                                                                                                      
                                                                                                                                              
    if (notnull == 1)                                                                                                                         
    {                                                                                                                                         
        inv <- 1/d                                                                                                                           
     } else {                                                                                                                                  
        inv <- solve(diag(d))                                                                                                                
     }                                                                                                                                         
                                                                                                                                              
    if (notnull != m)                                                                                                                         
    {                                                                                                                                         
        inv <- cbind(inv, matrix(0, nrow=notnull, ncol=(m - notnull)))                                                                       
        inv <- rbind(inv, matrix(0, nrow=(m-notnull), ncol=m))                                                                               
    }                                                                                                                                         
                                                                                                                                              
    mp <- s$v%*%inv%*%t(s$u)                                                                                                                  
                                                                                                                                              
    mp[abs(mp) < Eps] <- 0                                                                                                                   
    return(mp)                                                                                                                               
}                                                                                                                                            
                                                                                                                                              
                                                                                                                                              
                                                                                                                                              
                                                                                                                                              
                                                                                                                                              
                                                                                                                                              
distance2 <- function (x1, x2) {                                                                                                              
   temp <- x1 - x2                                                                                                                            
   sum(temp * temp)                                                                                                                           
}                                                                                                                                             
                                                                                                                                              
                                                                                                                                              
nnmf_als <-function(x, k, maxiter, eps)                                                                                      
{                                                                                                                                             
                                                                                                                                              
    print_iter <- 50 # iterations between print                                                                                               
                                                                                                                                              
    x <- as.matrix(x)                                                                                                                         
                                                                                                                                              
    if (any(!is.finite(x)))                                                                                                                   
                                                                                                                                              
        stop("infinite or missing values in 'x'")                                                                                             
                                                                                                                                              
    dx <- dim(x)                                                                                                                              
                                                                                                                                              
    D <- dx[1L]                                                                                                                               
                                                                                                                                              
    N <- dx[2L]                                                                                                                               
    Xscale = sum(x)                                                                                                                           
                                                                                                                                              
                                                                                                                                              
    if (!D || !N)                                                                                                                            
                                                                                                                                              
        stop("0 extent dimensions")                                                                                                          
                                                                                                                                              
                                                                                                                                              
    W <- matrix(abs(rnorm(D * k)), D, k)                                                                                                      
                                                                                                                                              
    H <- matrix(abs(rnorm(k * N)), k, N)                                                                                                      
                                                                                                                                              
                                                                                                                                              
                                                                                                                                              
    Rscale = sum(W %*% H)                                                                                                                     
                                                                                                                                              
    sqrnorm = sqrt(Rscale / Xscale)                                                                                                           
                                                                                                                                              
    H = H / sqrnorm                                                                                                                           
    W = W / sqrnorm                                                                                                                           
                                                                                                                                              
    Xr_old = W %*% H                                                                                                                          
                                                                                                                                              
                                                                                                                                              
    for (iter in 1:maxiter) {                                                                                                                 
                                                                                                                                              
       W = x %*% t(mpinv(H%*%t(H)) %*% H)                                                                                                     
       W = (W>0) * W                                                                                                                          
       W = W/(t(matrix(rep(colSums(W),D),ncol(W),nrow(W))) + eps)                                                                             
                                                                                                                                              
       H =t(W %*% mpinv(t(W) %*% W)) %*% x                                                                                                    
                                                                                                                                              
       H = H * (H>0)                                                                                                                          
                                                                                                                                              
       if (iter%% print_iter ==0) {                                                                                                           
                                                                                                                                              
           Xr = W %*% H                                                                                                                       
                                                                                                                                              
           diff = sum(abs(Xr_old-Xr))                                                                                                         
                                                                                                                                              
           Xr_old = Xr                                                                                                                        
                                                                                                                                              
           eucl_dist = distance2(x, W %*% H)                                                                                                  
                                                                                                                                              
           errorx = mean(abs(x - W %*% H)) / mean(x)                                                                                          
                                                                                                                                              
           cat('Iter = ', iter , "\t")                                                                                                        
                                                                                                                                              
           cat('relative error = ', errorx, "\t")                                                                                             
                                                                                                                                              
           cat('diff = ', diff, "\t")                                                                                                         
                                                                                                                                              
           cat('eucl dist = ', eucl_dist, "\n")                                                                                               
                                                                                                                                              
           if (errorx < 10e-6) {                                                                                                              
                                                                                                                                              
                  cat("Execution finishes at iteration = ", iter, "\n")                                                                       
                                                                                                                                              
                  break                                                                                                                       
                                                                                                                                              
           }                                                                                                                                  
                                                                                                                                              
        }                                                                                                                                     
                                                                                                                                              
     }                                                                                                                                        
     z <-c (list(W = W,H = H))                                                                                                                
                                                                                                                                              
     z                                                                                                                                        
                                                                                                                                              
}                                                                                                                                             
                      

nnmf_mm <-function(x, k, maxiter, eps)                                                                                    
{                                                                                                                                       
                                                                                                                                        
    print_iter <- 50 # iterations between print                                                                                        
                                                                                                                                        
    x <- as.matrix(x)                                                                                                                   
                                                                                                                                        
    if (any(!is.finite(x)))                                                                                                             
                                                                                                                                        
        stop("infinite or missing values in 'x'")                                                                                       
                                                                                                                                        
    dx <- dim(x)                                                                                                                        
                                                                                                                                        
    n <- dx[1L]                                                                                                                         
                                                                                                                                        
    m <- dx[2L]                                                                                                                         
                                                                                                                                        
    if (!n || !m)                                                                                                                       
                                                                                                                                        
        stop("0 extent dimensions")                                                                                                     
                                                                                                                                        
                                                                                                                                            
    W <- matrix(abs(rnorm(n * k)), n, k)                                                                                                
                                                                                                                                        
    H <- matrix(abs(rnorm(k * m)), k, m)                                                                                                
                                                                                                                                        
    Xr_old = W %*% H                                                                                                                    
                                                                                                                                        
    for (iter in 1:maxiter) {                                                                                                           
                                                                                                                                        
       H = H * (t(W) %*% x) / ((t(W) %*% W) %*% H + eps)                                                                                
                                                                                                                                        
       W = W * t(H %*% t(x)) / (W %*% (H %*% t(H)) + eps)                                                                               
                                                                                                                                        
       if (iter%% print_iter ==0) {                                                                                                     
                                                                                                                                        
           Xr = W %*% H                                                                                                                 
                                                                                                                                        
           diff = sum(abs(Xr_old-Xr))                                                                                                   
                                                                                                                                        
           Xr_old = Xr                                        
                                                     
           eucl_dist = distance2(x, W %*% H)
                    
           errorx = mean(abs(x - W %*% H)) / mean(x)                                                                                    
                                                                                                                                        
           cat('Iter = ', iter , "\t")                                                                                                 
                                                                                                                                        
           cat('relative error = ', errorx, "\t")    

           cat('diff = ', diff, "\t")                                                                                     

           cat('eucl dist = ', eucl_dist, "\n")                                                
                                                                                                                                        
           if (errorx < 10e-6) {                                                                                                        
                                                                                                                                        
                  cat("Execution finishes at iteration = ", iter, "\n")                                                                  
                                                                                                                                        
                  break                                                                                                                 
                                                                                                                                        
            }                                                                                                                            
                                                                                                                                        
         }

       }                                                                                                                                
                                                                                                                                        
     z <-c (list(W = W,H = H))                                                                 
     z                                                                                                                                  
}






 

                                                                                                                     
                                                                                                                                              
                                                                                                                                              
 

                                                                                                                                         
                                                                                                                                         
nnmf_prob <- function(x, k, maxiter, eps = 100*.Machine$double.eps)                                                                      
{                                                                                                                                        
                                                                                                                                         
  print_iter=50;                                                                                                                         
  powers=1.5+(2.5-1.5)*((1:maxiter)-1)/(maxiter-1);                                                                                      
                                                                                                                                         
  D=dim(x)[1L];                                                                                                                          
  N=dim(x)[2L];                                                                                                                          
                                                                                                                                         
  X_factor = sum(x);                                                                                                                     
  X_org = x;                                                                                                                             
  x=x/X_factor;                                                                                                                          
                                                                                                                                         
  W<-matrix(abs(rnorm(D*k)), D,k)                                                                                                        
                                                                                                                                         
  W = W / t(matrix(rep(colSums(W),D), ncol(W),nrow(W)))                                                                                  
                                                                                                                                         
  H<-matrix(abs(rnorm(k*N)),k,N)                                                                                                         
                                                                                                                                         
  H= H / (matrix(rep(rowSums(H),N), nrow(H),ncol(H)))                                                                                    
                                                                                                                                         
                                                                                                                                         
  P=matrix(rep(1),k,1)                                                                                                                   
                                                                                                                                         
  P=P/sum(P)                                                                                                                             
                                                                                                                                         
  W1=W                                                                                                                                   
  H1=H                                                                                                                                   
                                                                                                                                         
                                                                                                                                         
                                                                                                                                         
  Xr_old = W%*%H                                                                                                                         
  for (iter in 1:maxiter) {                                                                                                              
                                                                                                                                         
    Qnorm = (W %*% diag(diag(P),k,k)) %*%H                                                                                               
                                                                                                                                         
    for (j in 1:k) {                                                                                                                     
      Q<- ( t(t(W[,j])) %*% H[j,] * P[j] ) / (Qnorm + eps)                                                                               
      XQ = x * Q                                                                                                                         
                                                                                                                                         
                                                                                                                                         
      dummy =  rowSums(XQ)                                                                                                               
      W1[,j]=t(dummy / sum(dummy))                                                                                                       
                                                                                                                                         
      dummy = colSums(XQ)                                                                                                                
      H1[j,]=(dummy / sum(dummy))                                                                                                        
                                                                                                                                         
    }                                                                                                                                    
    W=W1                                                                                                                                 
    H=H1                                                                                                                                 
                                                                                                                                         
    if (iter%% print_iter ==0) {                                                                                                         
                                                                                                                                         
           Xr = W %*% H                                                                                                                  
                                                                                                                                         
           diff = sum(abs(Xr_old-Xr))                                                                                                    
                                                                                                                                         
           Xr_old = Xr                                                                                                                   
                                                                                                                                         
           eucl_dist = distance2(x, W %*% H)                                                                                             
                                                                                                                                         
           errorx = mean(abs(x - W %*% H)) / mean(x)                                                                                     
                                                                                                                                         
           cat('Iter = ', iter , "\t")                                                                                                   
                                                                                                                                         
           cat('relative error = ', errorx, "\t")                                                                                        
                                                                                                                                         
           cat('diff = ', diff, "\t")                                                                                                    
                                                                                                                                         
           cat('eucl dist = ', eucl_dist, "\n")                                                                                          
                                                                                                                                         
           if (errorx < 10e-6) {                                                                                                           
                                                                                                                                         
                  cat("Execution finishes at iteration = ", iter, "\n")                                                                  
                                                                                                                                         
                  break                                                                                                                  
                                                                                                                                         
            }                                                                                                                            
                                                                                                                                         
       }                                                                                                                                 
                                                                                                                                         
  }                                                                                                                                      
                                                                                                                                         
W =  W%*%diag(diag(sqrt(P)),k,k) *X_factor                                                                                               
H = diag(diag(sqrt(P)),k,k) %*%H                                                                                                         
                                                                                                                                         
z<- c(list(W=W, H=H))                                                                                                                    
                                                                                                                                         
z                                                                                                                                        
                                                                                                                                         
}                                                                                                                                        
                                                                                                                                         
                  


                                                                                                                        
                                                                                 
nnmf <- function(x, k, method = 'nnmf_mm', maxiter = 1000, eps=2.2204e-016)
{
      if (method == 'nnmf_als') {
            cat('Alternating Least Squares Algorithm', '\n')
            nnmf_als(x, k, maxiter, eps)
            

      }         

      else if (method == 'nnmf_prob') {
            cat('Multinomial Algorithm', '\n')
            nnmf_prob(x, k, maxiter, eps)
            

      }         


      else {
           cat('Multiplicative Update Algorithm', '\n')
           nnmf_mm(x, k, maxiter, eps)
      }

}