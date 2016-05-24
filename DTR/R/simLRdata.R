###################################################
### Reference:
### Kidwell KM, Wahed AS: Weighted log-rank statistic to compare shared-path  
### adaptive treatment strategies. Biostatistics. 14(2):299-312, 2013
###################################################

###################################################
### code chunk number 1: chunklibraries
###################################################

require(survival)

###################################################
### code chunk number 2: chunksimulatedata
###################################################

simLRdata<-function(n, # Number of subjects assigned to one arm
                           max.c, # Censoring time C is generated from Uniform(min=0,max=max.c)
                           pi.x, # First-stage indicator X is generated from Bernoulli(pi.x); X=0 for assignment to A1 and X=1 for assignment to A2
                           pi.r, # Remission/consent indicator R is generated from Bernoulli(pi.r); R=0 for nonresponders and R=1 for responders
                           pi.z, # Second-stage indicator Z is generated from Bernoulli(pi.z); Z=0 for assignment to B1 and Z=1 for assignment to B2
                           mean.NR.1, # For nonresponders (R=0) assigned to A1 at first stage, survival time T.NR.1 is drawn from exponential(1/mean.NR.1)
                           mean.NR.2, # For nonresponders (R=0) assigned to A2 at first stage, survival time T.NR.2 is drawn from exponential(1/mean.NR.2)
                           mean.R.1, # For responders (R=1) assigned to A1 at first stage, time to response T.R.1 is drawn from exponential(1/mean.R.1)
                           mean.R.2, # For responders (R=1) assigned to A2 at first stage, time to response T.R.2 is drawn from exponential(1/mean.R.2)
                           mean.RE.11, # For responders (R=1) assigned to A1 at first stage and B1 at second stage, time from response to event is drawn from exponential(1/mean.RE.11)
                           mean.RE.12, # For responders (R=1) assigned to A1 at first stage and B2 at second stage, time from response to event is drawn from exponential(1/mean.RE.12)
                           mean.RE.21, # For responders (R=1) assigned to A2 at first stage and B1 at second stage, time from response to event is drawn from exponential(1/mean.RE.21)
                           mean.RE.22 # For responders (R=1) assigned to A2 at first stage and B2 at second stage, time from response to event is drawn from exponential(1/mean.RE.22)
) { 

  #Function to test if a number is an integer or not.
  is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    
  #Function to check if a number is a valid probability value or not
  is.probability <- function(x) if(x>=0 && x<=1) TRUE else FALSE
    
  #Step 0: Check input errrors
  if (is.null(n)) stop("n can not be empty") 
  if (is.character(n)) stop("n has to be numeric")
  if (n<=0 || !is.wholenumber(n)) stop("n must be a positive integer")
  
  if (is.null(max.c)) stop("max.c can not be empty") 
  if (is.character(max.c)) stop("max.c has to be numeric")
  if (max.c<=0) stop("max.c must be a positive value") 
  
  if (is.null(pi.x)) stop("pi.x can not be empty") 
  if (is.character(pi.x)) stop("pi.r has to be numeric")
  if (!is.probability(pi.x)) stop("pi.x must be a valid probability value")
  
  if (is.null(pi.r)) stop("pi.r can not be empty") 
  if (is.character(pi.r)) stop("pi.r has to be numeric")
  if (!is.probability(pi.r)) stop("pi.r must be a valid probability value")
  
  if (is.null(pi.z)) stop("pi.z can not be empty") 
  if (is.character(pi.z)) stop("pi.z has to be numeric")  
  if (!is.probability(pi.z)) stop("pi.z must be a valid probability value")
  
  if (is.null(mean.NR.1)) stop("mean.NR.1 can not be empty") 
  if (is.character(mean.NR.1)) stop("mean.NR.1 has to be numeric")
  if (mean.NR.1<=0) stop("mean.NR.1 must be positive for the exponential distribution")
  
  if (is.null(mean.NR.2)) stop("mean.NR.2 can not be empty") 
  if (is.character(mean.NR.2)) stop("mean.NR.2 has to be numeric")
  if (mean.NR.2<=0) stop("mean.NR.2 must be positive for the exponential distribution")
  
  if (is.null(mean.R.1)) stop("mean.R.1 can not be empty") 
  if (is.character(mean.R.1)) stop("mean.R.1 has to be numeric")
  if (mean.R.1<=0) stop("mean.R.1 must be positive for the exponential distribution")
  
  if (is.null(mean.R.2)) stop("mean.R.2 can not be empty") 
  if (is.character(mean.R.2)) stop("mean.R.2 has to be numeric")
  if (mean.R.2<=0) stop("mean.R.2 must be positive for the exponential distribution") 
  
  if (is.null(mean.RE.11)) stop("mean.RE.11 can not be empty") 
  if (is.character(mean.RE.11)) stop("mean.RE.11 has to be numeric")
  if (mean.RE.11<=0) stop("mean.RE.11 must be positive for the exponential distribution")
  
  if (is.null(mean.RE.12)) stop("mean.RE.12 can not be empty") 
  if (is.character(mean.RE.12)) stop("mean.RE.12 has to be numeric")
  if (mean.RE.12<=0) stop("mean.RE.12 must be positive for the exponential distribution")
  
  if (is.null(mean.RE.21)) stop("mean.RE.21 can not be empty") 
  if (is.character(mean.RE.21)) stop("mean.RE.21 has to be numeric")
  if (mean.RE.21<=0) stop("mean.RE.21 must be positive for the exponential distribution")
  
  if (is.null(mean.RE.22)) stop("mean.RE.22 can not be empty") 
  if (is.character(mean.RE.22)) stop("mean.RE.22 has to be numeric")
  if (mean.RE.22<=0) stop("mean.RE.22 must be positive for the exponential distribution")   
  
  #Step 1
  #Generate censoring
  C <- runif(n,min=0,max=max.c)
  
  #Generate indicator for first-stage therapies A1/A2
  X <- rbinom(n, 1, pi.x) # X=0 for A1, 1 for A2
  
  #Generate Remission/Consent indicator
  R <- rbinom(n, 1, pi.r)
      
  #Generate B treatment indicator 
  #Only for R=1. Because Z is only supposed to be generated for responders
  #Non-responders will not participate in the second randomization 
  #For non-responders, Z is automatically set to 0
  Z <- rep(0,n)
  Z[which(R==1)] <- rbinom(length(which(R==1)),1,pi.z) # Z=0 for B1, 1 for B2
      
  #Step 2
  #When X=0 and R=0, a survival time T.NR.1 was drawn from exponential(1/mean.NR.1)
  T.NR.1 <- rep(0,n)
  T.NR.1[which(X==0 & R==0)] <- rexp(length(which(X==0 & R==0)),rate=1/mean.NR.1)
  #When X=1 and R=0, a survival time T.NR.2 was drawn from exponential(1/mean.NR.2)
  T.NR.2 <- rep(0,n)
  T.NR.2[which(X==1 & R==0)] <- rexp(length(which(X==1 & R==0)),rate=1/mean.NR.2)
      
  #When X=0 and R=1, time to response T.R.1 was drawn from exponential(1/mean.R.1)
  T.R.1 <- rep(0,n)
  T.R.1[which(X==0 & R==1)] <- rexp(length(which(X==0 & R==1)),rate=1/mean.R.1)
  #When X=1 and R=1, time to response T.R.2 was drawn from exponential(1/mean.R.2)
  T.R.2 <- rep(0,n)
  T.R.2[which(X==1 & R==1)] <- rexp(length(which(X==1 & R==1)),rate=1/mean.R.2)
  
  #Generate observed TR
  TR <- rep(0, n)
  TR <- T.R.1 + T.R.2
  
  #When X=0, R=1 and Z=0, time from response to event T.RE.11 was drawn from exponential(1/mean.RE.11)
  T.RE.11 <- rep(0,n)
  T.RE.11[which(X==0 & R==1 & Z==0)] <- rexp(length(which(X==0 & R==1 & Z==0)),rate=1/mean.RE.11)
  #When X=0, R=1 and Z=1, time from response to event T.RE.12 was drawn from exponential(1/mean.RE.12)
  T.RE.12 <- rep(0,n)
  T.RE.12[which(X==0 & R==1 & Z==1)] <- rexp(length(which(X==0 & R==1 & Z==1)),rate=1/mean.RE.12)
  #When X=1, R=1 and Z=0, time from response to event T.RE.21 was drawn from exponential(1/mean.RE.21)
  T.RE.21 <- rep(0,n)
  T.RE.21[which(X==1 & R==1 & Z==0)] <- rexp(length(which(X==1 & R==1 & Z==0)),rate=1/mean.RE.21)
  #When X=1, R=1 and Z=1, time from response to event T.RE.22 was drawn from exponential(1/mean.RE.22)
  T.RE.22 <- rep(0,n)
  T.RE.22[which(X==1 & R==1 & Z==1)] <- rexp(length(which(X==1 & R==1 & Z==1)),rate=1/mean.RE.22)
      
  #Step 3
  #Generate potential survival times T11_star, T12_star, T21_star, and T22_star
  T11_star <- T.R.1 + T.RE.11
  T12_star <- T.R.1 + T.RE.12
  T21_star <- T.R.2 + T.RE.21
  T22_star <- T.R.2 + T.RE.22
  
  #Step 4
  #Generate observed survival time T
  T <- (1-X)*((1-R)*T.NR.1 + R*(1-Z)*T11_star + R*Z*T12_star) + X*((1-R)*T.NR.2 + R*(1-Z)*T21_star + R*Z*T22_star)

  #Generate observed time U
  U <- pmin(T,C)
  
  #Generate censoring indicator
  #delta=0 for censoring and delta=1 for event
  delta <- as.numeric(T<C)
      
  #Step 5
  #Combine data
  data <- data.frame(X,TR,R,Z,U,delta)
  return(data)
  
}

