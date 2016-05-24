###################################################
### Reference:
### Guo X, Tsiatis AA: A weighted risk set estimator for survival distributions in two-stage 
### randomization designs with censored survival data. Int. J. Biostatistics 1:1-15, 2005
###################################################

###################################################
### code chunk number 1: chunklibraries
###################################################

require(survival)

###################################################
### code chunk number 2: chunksimulatedata
###################################################

simWRSEdata<-function(n, # Number of subjects assigned to one arm
                        max.c, # Censoring time C is generated from Uniform(min=0,max=max.c)
                        pi.r, # Remission/consent indicator R is generated from Bernoulli(pi.r); R=0 for nonresponders and R=1 for responders
                        pi.z, # Second-stage indicator Z is generated from Bernoulli(pi.z); Z=0 for assignment to B1 and Z=1 for assignment to B2
                        mean.T0, # For nonresponders (R=0), survival time T0 is drawn from exponential(1/mean.T0)
                        mean.TR, # For responders (R=1), time to response TR is drawn from exponential(1/mean.TR); for nonresponders (R=0), TR=0
                        mean.T1, # For assigned to B1 (Z=0), time T1_star is drawn from exponential(1/mean.T1)
                        mean.T2 # For assigned to B2 (Z=1), time T2_star is drawn from exponential(1/mean.T2)
) { 

  #Functions needed to check errors:
  #Function to test if a number is an integer or not.
  is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  #Function to check if a number is a valid probability value or not
  is.probability<-function(x) if(x>=0 && x<=1) TRUE else FALSE
  
  #Step 0: Check input errors
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
  
  if (is.null(mean.T0)) stop("mean.T0 can not be empty")
  if (is.character(mean.T0)) stop("mean.T0 has to be numeric")
  if (mean.T0<=0) stop("mean.T0 must be positive for the exponential distribution")
  
  if (is.null(mean.TR)) stop("mean.TR can not be empty")
  if (is.character(mean.TR)) stop("mean.TR has to be numeric")
  if (mean.TR<=0) stop("mean.TR must be positive for the exponential distribution")
  
  if (is.null(mean.T1)) stop("mean.T1 can not be empty")
  if (is.character(mean.T1)) stop("mean.T1 has to be numeric")
  if (mean.T1<=0) stop("mean.T1 must be positive for the exponential distribution")
  
  if (is.null(mean.T2)) stop("mean.T2 can not be empty")
  if (is.character(mean.T2)) stop("mean.T2 has to be numeric")    
  if (mean.T2<=0) stop("mean.T2 must be positive for the exponential distribution") 
  
  #Step 1
  #Generate censoring
  C<-runif(n,min=0,max=max.c)
  
  #Generate Remission/Consent indicator
  R<-rbinom(n,1,pi.r)
  
  #Generate B treatment indicator 
  #Only for R=1. Because Z is only supposed to be generated for responders
  #Non-responders will not participate in the second randomization 
  #For non-responders, Z is automatically set to 0
  Z<-rep(0,n)
  Z[which(R==1)]<-rbinom(length(which(R==1)),1,pi.z)
  
  #Step 2
  #When R=0, a survival time T0 was drawn from exponential(1/mean.T0)
  T0<-rep(0,n)
  T0[which(R==0)]<-rexp(length(which(R==0)),rate=1/mean.T0)
  
  #When R=1, time to response was drawn from exponential(1/mean.TR)
  TR<-rep(0,n)
  TR[which(R==1)]<-rexp(length(which(R==1)),rate=1/mean.TR)
  
  #For Z=0, a survival time T1_star was be generated from exponential(1/mean.T1)
  T1_star<-rep(0,n)
  T1_star[which(Z==0)]<-rexp(length(which(Z==0)),rate=1/mean.T1)
  
  #For Z=1, a survival time T2_star was be generated from exponential(1/mean.T2)   
  T2_star<-rep(0,n)
  T2_star[which(Z==1)]<-rexp(length(which(Z==1)),rate=1/mean.T2)
    
  #Step 3
  #Generate survival time T
  T<-(1-R)*T0+R*(TR+(1-Z)*T1_star+Z*T2_star)
  
  #Generate observed time U
  U<-pmin(T,C)
  
  #Generate censoring indicator
  #delta=0 for censoring and delta=1 for event
  delta<-as.numeric(T<C)

  #Step 4
  #Combine data
  data<-data.frame(TR,R,Z,U,delta)
  return(data)

}

