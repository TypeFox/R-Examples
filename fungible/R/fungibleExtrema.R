# Version: October 23, 2012
#                    R Code for fungibleExtrema
#
##  R function: fungibleExtrema
##
## Locate extrema of Fungible weights (i.e., find max (min) cos(b,a)
##  
##
##  Input Variables
##  R.X   p x p Predictor variable correlation matrix.
##  rxy   p x 1 Vector of predictor-criterion correlations.
##  r.yhata.yhatb  = correlation between alternate-weight (yhata) and 
##                   and  least squares (yhatb) composites.
##  Nstarts  Maximum number of (max) minimizations from random starting configurations.
##  MaxMin   Max = maximize cos(a,b); Min = minimize cos(a,b)
##
##  Output Variables
##  cos.ab  cosine between OLS and alternate weights
##  a     extrema of fungible weights
##  k     k weights
##  z     z weights
##  b     OLS weights
##  u     u weights
##  r.yhata.yhatb  correlation between yhat and ytilde
##  r.y.yhat       correlation between y and yhat
##  gradient   gradient at solution

        

fungibleExtrema <- function(R.X,rxy,r.yhata.yhatb,Nstarts=100,MaxMin="Max"){
 
  # auxiliary function definitions
  # vector norm
  norm <- function(x) x/as.numeric( sqrt(t(x) %*%x))
  # vector cosine 
  vec.cos <- function(x,y){ t(norm(x))%*%norm(y) }
  # vector length 
  lngth <- function(x) sqrt(t(x) %*%x)


##~~~~~~~~~~~~~Generate U matrix via Gram Schmidt~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
#    GenU.BAK <- function(mat,u){
#         ## 
#         p <- ncol(mat)
#         n <- nrow(mat)
#         oData <-matrix(0,n,p+1)
#         oData[,1]<-u
# 
#         for(i in 2:(p+1)){
#            oData[,i] <- resid(lm(mat[,(i-1)]~-1+oData[,1:(i-1)]))
#         }
# 
#         U<-oData[,2:(p+1)]
#         d <- diag(1/sqrt(diag(crossprod(U))))
#         U <- U%*%d
#         U
#     }#end GenU

##~~~~~~~~~~~~~Generate U matrix via QR Decomposition~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  GenU <- function(mat,u){
    qr.Q(qr(cbind(u,mat)))[,-1]
  }#end GenU

#~~~~~~~~~~~~~~~Compute Analytic Gradient~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

gradf <- function(sv){

   z<-sv[1:(NX-1)]       
   zpz <- t(z)%*%z 
   L <- sv[NX]        
   k <- r * u  + U %*% z * e     
   a <- V%*%Linv.sqrt%*%k
   apa <- t(a)%*%a
   ab <- as.numeric(t(a) %*%b.unit)
  
   dfdz <- (e * t(U) %*% Linv.sqrt %*%t(V)) %*% 
           (-as.numeric(ab*(apa)^-1.5) * a  + as.numeric((apa^-.5))*b.unit )
  
   
   if(MaxMin=="MIN") {
       dfdz <-   dfdz +4*L*(zpz-1)*z  
       dfdL<-    (zpz - 1)^2  
   }       
   if(MaxMin=="MAX") {
    dfdz <-   dfdz -4*L*(zpz-1)*z  
    dfdL<-   -(zpz - 1)^2    
   } 
    
   c(dfdz,dfdL)      
 }   

##~~~~~~~~~~~~Main Function~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  #sv = start values
  minv <- function(sv){
      
      z<-sv[1:(NX-1)]
      
      zpz <- t(z)%*%z
      L <- sv[NX]     
      
      k <- r * u + + U %*% z * sqrt(1-r^2)
            
      k.star <- Linv.sqrt%*%k
      a <- V %*% k.star  
      len.a<-lngth(a)
      f <-  (as.numeric(len.a)^-1)*t(a)%*%b.unit  + L*(zpz - 1)^2 
      
      if(MaxMin=="MAX") f<- (as.numeric(len.a)^-1)*t(a)%*%b.unit  - L*(zpz - 1)^2 
      if(MaxMin=="MIN") f<- (as.numeric(len.a)^-1)*t(a)%*%b.unit  + L*(zpz - 1)^2  
      f
   
   }  ## end minv    
     
 

  MaxMin <- toupper(MaxMin)  
  NX <- ncol(R.X)

  #OLS weights
  b <- crossprod(solve(R.X),rxy)
  b.unit <- norm(b)
  
  len.b <- sqrt(t(b)%*%b)
  r <- as.numeric(r.yhata.yhatb)
  e <- sqrt(1-r^2)

  VLV <- eigen(R.X)
  V <- VLV$vectors
  L <- diag(VLV$values)
  Linv<-solve(L)
  Linv.sqrt <- solve( sqrt(L))
 
  u <- (sqrt(L)%*%t(V)%*%b)/ as.numeric(sqrt(t(b)%*%R.X%*%b))

  mat <- matrix(rnorm(NX*(NX-1)),NX,NX-1)
  U <- GenU(mat,u)
    
     FNSCALE <- 1
     if(MaxMin=="MAX")  FNSCALE <- FNSCALE * -1  
    
     # declare minf /maxf
     minf <- 99*FNSCALE
     maxf <- -1    
     breakflag<-0
     iter <- 0
  
    while(iter < Nstarts){
      # start values
      sv <- c(rnorm(NX-1)/sqrt(NX-1),1000)
      
      tmp <- try(optim(par=sv,
                       fn=minv,
                       gr=gradf,
                       method="BFGS",
                       control=list(fnscale=FNSCALE,
                                    maxit=500,
                                    parscale=c(rep(1,NX-1),1))))
                
                
                 if(abs(tmp$value)>1) tmp$convergence<-1
                 if((FNSCALE*tmp$value <= FNSCALE*minf) & (tmp$convergence==0))
                 { 
                    iter <- iter + 1
                    fdelta <-minf-tmp$value
                    minf<-tmp$value
                    out<-tmp           
                    z<-out$par[1:(NX-1)]
                    k <- r * u + U %*% z * sqrt(1-r^2)
                    a <- a.tilde <- V %*% Linv.sqrt %*% k 
                                       
                    scaling.weight <- (t(rxy) %*% a)/(t(a)%*%R.X%*%a)               
                    a <- as.numeric(scaling.weight) * a  
                    
                    cat(c(iter," Current fnc val",minf,"\n"))
                    print(max(abs(gradf(out$par))))
                    if(max(abs(gradf(out$par)))<=1e-5) breakflag <- 1           
                 } # End tmp$minimum < minf
               
       if(breakflag) break
   }# End NLoop    
      

      if(sign(a[1])!=sign(a.tilde[1]))  a <- a*-1

      r.y.yhat <- sqrt( (t(b) %*% R.X %*%b) )
      
      cat("\n\n")
      if(MaxMin=="MAX") cat("\n\n  Maximizing cos(a,b):\n")
      if(MaxMin=="MIN") cat("\n\n  Minimizing cos(a,b):\n")        
      cat("  r(yhat.a,yhat.b) = ",round(r.yhata.yhatb,3),"\n")
	  #compute vector cosine 
      s<-vcos(a,b) 
      cat("  cos(a,b) = ",round(s,3),"\n")
      cat("  RSQb = ",round(r.y.yhat^2,3),"\n")
      cat("  RSQa = ",round((r.yhata.yhatb * r.y.yhat)^2,3),"\n")
      cat("  Relative loss = RSQb - RSQa = ",
             round( r.y.yhat^2-(r.yhata.yhatb * r.y.yhat)^2 ,3),"\n\n")
         
      ba.mat <- cbind(b,a)
      colnames(ba.mat) <- c("b","a")
      print(ba.mat)
      cat("\n\n")
     
    
      cat("Analytic Gradient at Solution\n")
      solution.gradient <- gradf(out$par)
      print(matrix(solution.gradient,NX,1))
      cat("\n")
         

         
   list(cos.ab = s,
        a = a ,
        k = k,
        z=z,
        b = b,
        u = u,
        r.yhata.yhatb=r.yhata.yhatb,
        r.y.yhatb=r.y.yhat,
        gradient=solution.gradient)
 } ## End of Fungible Function





# ##----------------------------------------##
# ##   EXAMPLE
# ##  This is Koopmnan's Table 2 Example
# ##----------------------------------------##
# 
# R.X <- matrix(c(1.00,  .69,  .49,  .39,
#                  .69, 1.00,  .38,  .19,
#                  .49,  .38, 1.00,  .27,
#                  .39,  .19,  .27, 1.00),4,4)
#                  
#                  
# b <- c(.39, .22, .02, .43)
# rxy <- R.X%*%b
# 
# OLSRSQ <- t(b)%*%R.X%*%b
# 
# #theta <- .02
# #r.yhata.yhatb <- sqrt( 1 - (theta)/OLSRSQ)
# 
# r.yhata.yhatb  <- .90
# output <- FungibleExtrema1(R.X, rxy, r.yhata.yhatb, Nstarts=500, MaxMin="Min")
# 
#  ## Scale to replicate Koopman
#       a<-output$a
#       a.old<-a
#       aRa<-t(a)%*%R.X%*%a
#       ##Scale a such that a' R a = .68659
#       ##vc =variance of composite
#       vc <- aRa
#       ##sf = scale factor
#       sf <- .68659/vc
#       a <- as.numeric(sqrt(sf))*a
#       cat("\nKoopman Scaling\n")
#       print(round(a,2))
# 
  
 


