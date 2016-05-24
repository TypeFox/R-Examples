###################################################
### Reference:
### Lunceford JK, Davidian M, Tsiatis AA: Estimation of survival distributions of treatment 
### policies in two-stage randomization designs in clinical trials. Biometrics 58:48-57, 2002
###################################################

###################################################
### code chunk number 1: chunklibraries
###################################################

#Libraries required
require(survival)

###################################################
### code chunk number 2: chunksimulatedata
###################################################

simLDTdata<-function(n, # Number of subjects assigned to one arm
                       max.c, # Censoring time C is generated from Uniform(min=0,max=max.c)
                       pi.r, # Remission/consent indicator R is generated from Bernoulli(pi.r); R=0 for nonresponders and R=1 for responders
                       pi.z, # Second-stage indicator Z is generated from Bernoulli(pi.z); Z=0 for assignment to B1 and Z=1 for assignment to B2
                       lambda, # For nonresponders (R=0), survival time T.star.lambda is drawn from exponential(lambda)
                       alpha, # For responders (R=1), T.star.alpha is drawn from exponential(alpha)
                       beta1, # For responders (R=1), T.star.11 is drawn from exponential(exp(beta1))
                       beta2, # For responders (R=1), T.star.12 is drawn from exponential(exp(beta1+beta2*T.star.11))
                       L=.Machine$double.xmax # Optional, restricted survival time L
) { 
  
  #Functions needed to check errors:
  #Function to test if a number is an integer or not.
  is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    
  #Function to check if a number is a valid probability value or not
  is.probability <- function(x) if(x>=0 && x<=1) TRUE else FALSE
    
  #Step 0: Check input errors:
  if (is.null(n)) stop("n can not be empty")
  if (is.character(n)) stop("n has to be numeric")
  if (n<=0 || !is.wholenumber(n)) stop("n must be a positive integer")
  
  if (is.null(max.c)) stop("max.c can not be empty")
  if (is.character(max.c)) stop("max.c has to be numeric")  
  if (max.c<=0) stop("max.c must be a positive value")
  
  if (is.null(pi.r)) stop("pi.r can not be empty")
  if (is.character(pi.r)) stop("pi.r has to be numeric")
  if (!is.probability(pi.r)) stop("pi.r must be a valid probability value")
  
  if (is.null(pi.z)) stop("pi.z can not be empty")
  if (is.character(pi.z)) stop("pi.z has to be numeric")     
  if (!is.probability(pi.z)) stop("pi.z must be a valid probability value")
  
  if (is.null(lambda)) stop("lambda can not be empty")
  if (is.character(lambda)) stop("lambda has to be numeric")
  if (lambda<=0) stop("Lambda must be positive for the exponential distribution")
  
  if (is.null(alpha)) stop("alpha can not be empty")
  if (is.character(alpha)) stop("alpha has to be numeric")
  if (alpha<=0) stop("alpha must be positive for the exponential distribution") 
  
  if (is.null(beta1)) stop("beta1 can not be empty")
  if (is.character(beta1)) stop("beta1 has to be numeric")
  
  if (is.null(beta2)) stop("beta2 can not be empty")
  if (is.character(beta2)) stop("beta2 has to be numeric") 
  
  if (is.null(L)) stop("L can not be empty")
  if (is.character(L)) stop("L has to be numeric")
  if (L<=0) stop("L must be a positive value")
        
  #Step 1
  #Generate censoring
  C <- runif(n,min=0,max=max.c)
  
  #Generate Remission/Consent indicator
  R <- rbinom(n,1,pi.r)
      
  #Generate B treatment indicator 
  #Only for R=1. Because Z is only supposed to be generated for responders
  #Non-responders will not participate in the second randomization 
  #For non-responders, Z is automatically set to 0
  #Z=0 for assignment to B1 and Z=1 for assignment to B2
  Z <- rep(0,n)
  Z[which(R==1)] <- rbinom(length(which(R==1)),1,pi.z)
         
  #Step 2
  #When R=0, a survival time T.star.lambda was drawn from exponential(lambda)
  T.star.lambda <- rep(0,n)
  T.star.lambda[which(R==0)] <- rexp(length(which(R==0)),rate=lambda)
      
  #When R=1, T.star.alpha was drawn from exponential(alpha)
  T.star.alpha <- rep(0,n)
  T.star.alpha[which(R==1)] <- rexp(length(which(R==1)),rate=alpha)
    
  #Generate T.star.11
  T.star.11 <- rep(0,n)
  T.star.11[which(R==1)] <- rexp(length(which(R==1)),rate=exp(beta1))
      
  #Generate T.star.12 
  T.star.12 <- rep(0,n)
  T.star.12[which(R==1)] <- apply(as.array(T.star.11[which(R==1)]), 1, function(x) rexp(1,rate=exp(beta1+beta2*x)))
      
  #Step 3
  #Generate T.11 and T.12 based on T.star.lambda, T.star.alpha, T.star.11, and T.star.12
  T.11 <- pmin((1-R)*T.star.lambda+R*(T.star.alpha+T.star.11),rep(L,n))
  T.12 <- pmin((1-R)*T.star.lambda+R*(T.star.alpha+T.star.12),rep(L,n))
    
  #Step 4
  #Generate survival time T
  T <- (1-R)*T.11+R*(1-Z)*T.11+R*Z*T.12
  
  #Generate observed time U
  U <- pmin(T,C)
  
  #Generate censoring indicator
  #delta=0 for censoring and delta=1 for event
  delta <- as.numeric(T<C)
      
  #Step 5
  #Combine data
  data <- data.frame(R,Z,U,delta)
  return(data)
  
}



