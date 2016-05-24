####functions for EM
####5-22-07
################################EM
### use the function "ns" in library(splines)

##formatting data : create data object suitable for implementation of EM 
Format.EM <- function(data.list,n,nmax,grid){
   data.f<-Format.data(data.list,n,nmax)
   Obs<-data.f[[1]]
   T<-data.f[[2]]
   N<-data.f[[3]]
   data<-TranMtoV(Obs,T,N)
   y<-data[,2]
   t<-data[,3]

   
   timeindex = floor(t*length(grid))+1
   result = list(y=y,curve=data[,1],n=n,timeindex=timeindex)
   return(result)
}

###EM
EM<-function(data.list,n,nmax,grids,M.EM,iter.num,r.EM,basis.EM,sig.EM){
 
  R.inv<-R.inverse(M.EM,grids)
## (a) formatting data for use in the EM routine
 data.obj <- Format.EM(data.list,n,nmax,grids) 
## (b) EM estimation of eigenfunctions (theta) and PC scores (alpha)
 EMest <- fpcaEM(data.obj,k=r.EM,df=M.EM, grid = grids, maxit = iter.num, tol = 0.001, pert = 0.01, sigma = sig.EM^2,basis.EM,R.inv)

 if(basis.EM=="ns"){
 
 B.basis <- cbind(1, ns(grids, df = M.EM)) ### spline basis (with constant) used in EM
 B.orthobasis <- qr.Q(qr(B.basis)) ### orthonormalized spline basis used in EM
 }

if(basis.EM=="poly"){
     
    lmin<-0
    lmax<-1 
    delta<-(lmax-lmin)/(M.EM-3)
    knots<-c(seq(lmin,lmax,by=delta),lmax+delta,lmax+2*delta)
    bs<-apply(matrix(grids),1, BC.orth,knots=knots,R.inv=R.inv)
    B.orthobasis<-t(bs)*sqrt(grids[2]-grids[1])
} 

EMeigenmat <- B.orthobasis %*% EMest$theta ### (unscaled) eigenfunctions in original time scale
EMeigen.est <- EMeigenmat%*% diag(1/sqrt(diag(t(EMeigenmat)%*%EMeigenmat))) ###  normalization

## (c) Crude ''estimated'' of eigenvalues from normalization
EMDinv <- EMest$alpha[[3]]  ## estimate of D^{-1} from EM algorithm
##EMgamma <- EMeigenmat %*% solve(EMDinv) %*% t(EMeigenmat)
EMgamma <- EMeigenmat %*% solve(EMDinv, t(EMeigenmat))
### estimate of Gamma

####EMgamma.svd <- svd(EMgamma)  ## svd of Gamma (* length(grid))
####EMeigenvec.est <- EMgamma.svd$u[,1:r.EM]
### first r.EM eigenvectors of Gamma
#####EMeigenval.est <- EMgamma.svd$d[1:r.EM]/length(grids)


EMgamma.svd <- eigen(EMgamma,symmetric=TRUE)     ###use eigen with symmetric=TRUE
### first r.EM eigenvectors of Gamma
EMeigenval.est <- EMgamma.svd$values[1:r.EM]/length(grids)
### first r.EM eigenvectors of Gamma
EMeigenvec.est <- EMgamma.svd$vectors[,1:r.EM]
### estimated eigenvalues
EMsigma.est = sqrt(EMest$sigma)   ##estimated sigma

##result
result<-list(EMeigenvec.est,EMeigenval.est,EMsigma.est)
return(result)
 
}



### Code (sent by Gareth James) for implementation of EM algorithm in R
### Uses the procedure of James, Hastie and Sugar (2000)
### Uses the ``splines'' package and functions`bs' and 'ns' for representation of the 
### eigenfunctions

library(splines)

calclike <- function(y, sigma, Dinv, theta, B, curve){  

### calculates the log likelihood of data	

      like <- 0	
      N <- length(table(curve))	
      for(i in 1:N){       
           X <- B[curve == i, ] %*% theta		
           C.my <- X %*% solve(Dinv, t(X)) + sigma * diag(dim(X)[1])
           C.my<-(C.my+t(C.my))/2 
      like <- like - t(y[curve == i]) %*% solve(C.my, y[curve == i])/2 - sum(logb(eigen(C.my,symmetric=TRUE, only.values=TRUE)$value))/2	
      }
#     print(paste("Like =", like) 
      return(like)
}

fpcaEM <- function(obj, k = 2, df = 5, grid = seq(0.01, 1, length = 100), maxit = 50, tol = 0.001, pert = 0.01, sigma = 1, basis.method="ns",R.inv=NULL){

## computes the MLEs of alpha (PC score), theta (eigenfunctions represented in spline basis)

     y <- obj$y                 
### data represented as a single vector (by stacking the observations for different curves)

     timeindex <- obj$timeindex             
### a vector with the index of location (on the grid) of each point in y 

     curve <- obj$curve                   
### a vector with the curve number corresponding to each value of y 

     N <- obj$n                           
### number of observations (= n, in our notation)

 if(basis.method=="ns"){
   #  print("ns")
     B <- cbind(1, ns(grid, df = df))      
### spline basis evaluated at the grid locations, 
#### df = d.f. for the splines (= M in our notation); grid = a vector of possible time points


     B.orth <- qr.Q(qr(B))  #### orthonormalizing the columns of B
     B <- B.orth[timeindex, ]  ### evaluating the basis at the observed times
 } 
 
 if(basis.method=="poly"){
  #  print("poly")
 ###R.inv<-R.inverse(df)
    lmin<-0
    lmax<-1 
    delta<-(lmax-lmin)/(df-3)
    knots<-c(seq(lmin,lmax,by=delta),lmax+delta,lmax+2*delta)
    bs<-apply(matrix(grid),1, BC.orth,knots=knots,R.inv=R.inv)
    bs<-t(bs)*sqrt(grid[2]-grid[1])
    B<-bs[timeindex,]
    }

     R.old <- 0    
     R.new <- 1   
     ind <- 1  
     theta.zero <- solve(t(B) %*% B, t(B)) %*% y 
     y <- y - B %*% theta.zero  

     alpha <- list(alpha = matrix(1, N, k), alphaprod = array(1, c(N, k, k)))   
### k = dimension of alpha (= r, in our notation)

     theta <- as.matrix(init(y, timeindex, curve, B, k, pert)$theta)  
     alpha <- getemalpha(y, curve, theta, B, alpha)  
     while(abs(R.old - R.new)/R.new > tol & (ind < maxit)){    
          sigma <- getsigma(y, curve, B, theta, alpha, sigma)   
          theta <- getemtheta(y, curve, alpha, B, theta) 
          alpha <- getemalpha(y, curve, theta, B, alpha, sigma = sigma) 
          R.old <- R.new 
          R.new <- sum((y - (B %*% theta * alpha$alpha[curve, ]) %*% rep(1, k))^2)
#         print(R.new) 
          R.new
          calclike(y, sigma, alpha$Dinv, theta, B, curve) 
          ind <- ind + 1
     }  
     temp <- svd(theta)  
     result = list(alpha = alpha, theta = theta, B = B, theta.zero = theta.zero, sigma = sigma)
     return(result)
}

getemalpha <-function(y, curve, theta, B, alphaobj, sigma = 1){ 

### computes EM update for alpha (E-step of EM)

     alpha <- alphaobj$alpha  
     alphaprod <- alphaobj$alphaprod
     n <- length(table(curve))
     N <- dim(alpha)[1]  
     if(dim(alpha)[2] > 1) {
         alphasq <- apply(alphaprod, 1, diag)
         Dinv <- N * diag(as.vector(1/(alphasq %*% rep(1, N))))  
     }  
     else 
         Dinv <- as.matrix(1/mean(alphaprod))  
         for(i in 1:n) {  
               X <- B[curve == i, ] %*% theta 
               Calpha <- solve(sigma * Dinv + t(X) %*% X)  

               alpha[i, ] <- Calpha %*% t(X) %*% y[curve == i]              
###   hat{alpha_i} (eq. (27))

               alphaprod[i,  ,  ] <- sigma * Calpha + alpha[i,  ] %*% t(alpha[i,  ])      
### hat{alpha_i alapha_i^T}  (eq. (28))
     }  
     result = list(alpha = alpha, alphaprod = alphaprod, Dinv = Dinv)
     return(result)
}

getemtheta <- function(y, curve, alphaobj, B, theta, tol = 0.0001){ 

### computes EM update for theta (M step)

      q <- dim(B)[2]
      R.old <- 1
      R.new <- 0
      alpha <- alphaobj$alpha
      alphaprod <- alphaobj$alphaprod

      k <- dim(alpha)[2]              
### = r, in our notation

      while(abs(R.old - R.new)/R.new > tol) {	
           for(j in 1:k) {
               ind <- rep(1, k)   
               ind[j] <- 0       
               tempy <- alpha[curve, j] * y - ((B %*% theta) * alphaprod[curve,  ,  j]) %*% ind   
               tempX <- B * sqrt(alphaprod[curve, j, j])      
               theta[, j] <- solve(t(tempX) %*% tempX, t(B)) %*% tempy		   
           }    
           R.old <- R.new  
           R.new <- sum((y - ((B %*% theta) * alpha[curve,  ]) %*% rep(1, k))^2)
       }	    
       return(theta)
}

getsigma <- function(y, curve, B, theta, alpha, sigma){	

### EM update for sigma^2 (M step)

      tempsigma <- 0  
      Dinv <- alpha$Dinv   
      N <- dim(alpha$alpha)[1]
      fit <- NULL
      for(i in 1:N){   
          X <- B[curve == i,  ] %*% theta  
          Calpha <- solve(Dinv + t(X) %*% X/sigma)  
          fit <- c(fit, B[curve == i,  ] %*% theta %*% alpha$alpha[i,])
          tempsigma <- tempsigma + sum(diag(X %*% Calpha %*% t(X)))
      }	
      sigma <- (sum((y - fit)^2) + tempsigma)/length(y)
#     print(paste("sigma =", sigma))	
      return(sigma)
}

init <- function(y, timeindex, curve, B, k, pert = 0){     

### initial values of theta and gamma

      tab <- table(curve)  
      N <- length(tab) 
      s <- c(0, cumsum(tab)) 
      q <- dim(B)[2]   
      gamma <- matrix(0, N, q) 
      for(i in 1:N) {
         X <- B[(s[i] + 1):s[i + 1],  ]  
         tempy <- y[(s[i] + 1):s[i + 1]]     
         gamma[i,  ] <- solve(t(X) %*% X + pert * diag(q), t(X)) %*% tempy
      }	
      theta <- prcomp(gamma)$rotation[, 1:k]
      result = list(theta = theta, gamma = gamma)
      return(result)
}