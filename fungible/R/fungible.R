##################################################################
##
##                    R Code for Fungible Weights
##
##  R function: fungible
##  Niels Waller
##  March 11, 2008
##  modified August 12, 2015
##
##  Input Variables
##  R.X   p x p Predictor variable correlation matrix.
##  rxy   p x 1 Vector of predictor-criterion correlations.
##  r.yhata.yhatb  = correlation between least squares (yhatb)
##                   and alternate-weight (yhata) composites.
##  sets  Number of returned sets of fungible weights.
##  print    logical, if TRUE then print 5-point summaries
##           of alternative weights
##
##  Output Variables
##  a     sets x p matrix of fungible weights
##  k     sets x p matrix of k weights
##  b     p x 1 vector of LS weights
##  u     p x 1 vector of u weights
##  r.yhata.yhatb  correlation between yhata and yhatb
##  r.y.yhatb       correlation between y and yhatb
##  cov.a   Expected covariance matrix for a
##  cor.a   Expected correlation matrix for a
####################################################################

fungible <- function(R.X,rxy,r.yhata.yhatb,sets,print=TRUE){


##~~~~~~~~~~~~~Generate U matrix via QR Decomposition~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   GenU <- function(mat,u){
     ## Generate U matrix via QR
     qr.Q(qr(cbind(u,mat)))[,-1]
   }#end GenU
   
   
   
  NX <- ncol(R.X)
  a.matrix <- k.matrix <- matrix(0,sets,NX)

  #OLS weights
  b <- crossprod(solve(R.X),rxy)
  r <- as.numeric(r.yhata.yhatb)

  VLV <- eigen(R.X)
  V <- VLV$vectors
  L <- diag(VLV$values)
  Linv.sqrt <- solve( sqrt(L))
  u.star <- t(V)%*%b

  u.circle <- sqrt(L) %*% u.star
  u <- u.circle/ as.numeric(sqrt((t(u.circle) %*% u.circle)))
 
  r.y.yhatb <- sqrt( (t(b) %*% R.X %*%b) )
   
  mat <- matrix(rnorm(NX*(NX-1)),NX,NX-1)
  U <- GenU(mat,u)   
 
  for(i in 1:sets){

      z <- rnorm((NX-1))
      z <- z / as.numeric( sqrt( t(z) %*% z))
               
      k <- r * u + + U %*% z * sqrt(1-r^2)
      k.star <- Linv.sqrt%*%k
      a <- V %*% k.star
      
#  scale a to minimize SSE_a
      s <- (t(rxy) %*% a)/(t(a)%*%R.X%*%a)
      
      a <- as.numeric(s) * a

      if(i==1) {
         cat("\n\nGenerating alternate weights . . . \n")
         r.yhata.yhatb <- (t(a) %*% R.X %*%b)/( sqrt((t(a)%*%R.X%*%a))*
                          sqrt(t(b)%*%R.X%*%b))
      }
      
      a.matrix[i,] <-a
      k.matrix[i,] <- k
   }
   cat("\n\n")
   cat("  r.yhata.yhatb = ",r.yhata.yhatb,"\n")
   cat("  RSQb = ",round(r.y.yhatb^2,3),"\n")
   cat("  RSQa = ",round((r.yhata.yhatb * r.y.yhatb)^2,3),"\n")
   cat("  Relative loss = RSQb - RSQa = ",
   round( r.y.yhatb^2-(r.yhata.yhatb * r.y.yhatb)^2 ,3),"\n")
   cat("  OLS b = ",t(round(b,3)),"\n\n")
   cat("\n")
   colnames(a.matrix) <- paste("a",1:NX,sep="")
   
   if(print){
    cat("\nSummary of generated alternate weights\n")
    print( apply(a.matrix,2,summary) )
    cat("\n")
   } 
   
    
# Compute Expected Moments
    G <- V%*%Linv.sqrt%*%U
    esq <- (1-r^2)
# Expected a 
    mn.a <-  r^2 * b
    cat("\n Expected a \n")
    print(mn.a)
    Ezsq <- 1/(NX-1) 

# Expected covariance matrix 
    cov.a <- as.numeric(r^2 * r.y.yhatb^2)* Ezsq * esq * G %*%t(G)            
             
                                   
    cat("\nExpected Covariance Matrix \n")
    print(cov.a)
    Dmat <- diag(1/sqrt(diag(cov.a)))
    cor.a <- Dmat %*% cov.a %*% Dmat
    cat("\nExpected Correlation Matrix \n")
    print(cor.a)
    
       
   list(a = a.matrix ,
        k = k.matrix,
        b = b,
        u = u,
        r.yhata.yhatb=r.yhata.yhatb,
        r.y.yhatb=r.y.yhatb,
        cov.a=cov.a,
        cor.a=cor.a)
 } ## End of Fungible 
############################################################

# 
# 
# 
# ##-------------------##
# ##   EXAMPLE
# ##  GRE/GPA Data
# ##-------------------##
# 
# R.X <- matrix(c(1.00,  .56,  .77,
#                  .56, 1.00,  .73,
#                  .77,  .73, 1.00),3,3)
# 
# rxy <- c(.39, .34, .38)
# 
# 
# b <- solve(R.X)%*%rxy
# theta <- .01
# OLSRSQ <- t(b)%*%R.X%*%b                  
# 
# r.yhata.yhatb <- sqrt( 1 - (theta)/OLSRSQ)
# 
# Nsets <- 50
# 
# 
# output <- fungible(R.X, rxy, r.yhata.yhatb, sets=Nsets, print=TRUE)
# 

    
