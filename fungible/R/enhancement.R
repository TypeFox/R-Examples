#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
# Author: Niels G Waller
# updated September 12, 2013
# updated March 11, 2011
# For a given R, generate b such that 
# b'r = c1, r'r = c2 
# Arguments:
# R = predictor correlation matrix
# br = model R-squared = b' r
# rr = sum of squared predictor/criterion corrs
# Values:
# b = vector of standardized regression coefficients
# r = vector of predictor/criterion correlations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
enhancement<- function(R,br, rr)
{ 
   Nvar <- nrow(R)
   lambda <- eigen(R)$values 
   V <- eigen(R)$vectors
   L <- matrix(rbind(lambda,lambda^2),2,Nvar)

   # is r'r in valid range
  lambda.p<-min(lambda)   
   if(rr < br*lambda.p){ 
     warning("\n\nFATAL ERROR \n*** rr must be >= ", br*lambda.p,"  ***\n") 
     return(NA)
   } 	
   
   # find b when enhancement at max value 
   atMax <- FALSE 
   
   if( identical(round(br - rr,12),round( br*(1-lambda.p),12))){
    vp <- V[,ncol(R)] #last eigenvector 
    b <- as.numeric(sqrt(br/lambda.p))*vp 
    r <- R %*% b 
    atMax <- TRUE
   }
   
   # find b when enhancement  not at max value 
   if(!atMax){
     # vectors b'r = c1, r'r = c2 c1 >=c2 
     g <- matrix(c(br,rr),2,1) 

     # compute the g-inverse of L 
     L.ginv <- t(L) %*% solve(L %*%t(L)) 
 
     # initialize delta.sq
     delta.sq <- rep(-9,Nvar) 
     # all elements of delta.sq must be >= 0 
     k<-1
     L.ginvg<- L.ginv %*% g
     diagNvarMinusLginvL <- (diag(Nvar) - L.ginv %*% L) 
     Nvar5 <- 1/(5*Nvar)
     while( min(delta.sq) < 0 ){
         h <- matrix(rnorm(Nvar,0, Nvar5),Nvar,1) 
         delta.sq <- L.ginvg  + diagNvarMinusLginvL %*% h
      } 
 
     signMatrix <- diag(sample(c(-1,1),Nvar,replace=TRUE))
     delta <- signMatrix%*%sqrt(delta.sq) 
     b <- V %*% delta 
  } ## end!atMax
  
   list(b=b, r=R%*%b)  
} ##### End of enhancement

#=========== Example 1 ===========#
# rm(output,r,b)
# R <- matrix(c( 1, .5, .25,
#                .5, 1,  .30,
#                .25, .30, 1), 3, 3) 
# 
# Rsq = .60
# 
# output<-enhancement(R,br=Rsq,rr=.40) 
# output
# 
# r<-output$r
# b<-output$b
#  
# #Standardized regression coefficients
# t(b) 
# #Predictor-criterion correlations"
# t(r) 
# 
# #Coefficient of determinations (b'r)"
# t(b) %*% r
# #Sum of squared correlations (r'r)"
# t(r) %*% r

#=========== Example 2 ===========#
#Nvar<- 100
#RSQ <- .6
#R <- matrix( .1, 100, 100); diag(R) <- 1
#rr <- mean(c(.6,(eigen(R)$val[Nvar] * .6)))
#cat("\nInput Rsq = ", RSQ, "\nInput r'r = ", rr, "\n")
#        
#output<-enhancement(R,br=RSQ,rr) 
#
#
#r<-output$r
#b<-output$b
#
#cat("\nReproduced Rsq = ", t(b) %*% r)
#cat("\nReproduced r'r = ", t(r) %*% r)
#
#cat("\nNumber of predictors = ", Nvar, "\n")


