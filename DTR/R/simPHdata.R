###################################################
### Reference:
### Tang X, Wahed AS: Comparison of treatment regimes with adjustment for auxiliary
### variables. Journal of Applied Statistics 38(12):2925-2938, 2011
###################################################

###################################################
### code chunk number 1: chunklibraries
###################################################

require(survival)

###################################################
### code chunk number 2: chunksimulatedata
###################################################

simPHdata<-function(n, # Number of subjects in the clinical trial
                       max.c, # Censoring time C is generated from Uniform(min=0,max=max.c)
                       pi.x, # First-stage indicator X is generated from Bernoulli(pi.x); X=0 for assignment to A1 and X=1 for assignment to A2
                       pi.z, # Second-stage indicator Z is generated from Bernoulli(pi.z); Z=0 for assignment to B1 and Z=1 for assignment to B2
                       lambda, # Baseline hazard
                       alpha, # Time to response TR is generated from exponential(alpha)
                       beta1, # Coefficient for first-stage indicator X
                       beta2, # Coefficient for R(t)
                       beta3, # Coefficient for XR(t)
                       beta4, # Coefficient for ZR(t)
                       beta5, # Coefficient for XZR(t)
                       gamma # Coefficient vector for covariates
) { 

  #Functions needed to check errors:
  #Function to test if a number is an integer or not.
  is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  #Function to check if a number is a valid probability value or not
  is.probability<-function(x) if(x>=0 && x<=1) TRUE else FALSE
  
  #Check input errors
  if (is.null(n)) stop("n can not be empty")
  if (is.character(n)) stop("n has to be numeric")
  if (n<=0 || !is.wholenumber(n)) stop("n must be a positive integer")
  
  if (is.null(max.c)) stop("max.c can not be empty")
  if (is.character(max.c)) stop("max.c has to be numeric")
  if (max.c<=0) stop("max.c must be a positive value")
  
  if (is.null(pi.x)) stop("pi.x can not be empty")
  if (is.character(pi.x)) stop("pi.x has to be numeric")
  if (!is.probability(pi.x)) stop("pi.x must be a valid probability value") 
  
  if (is.null(pi.z)) stop("pi.z can not be empty")
  if (is.character(pi.z)) stop("pi.z has to be numeric")
  if (!is.probability(pi.z)) stop("pi.z must be a valid probability value") 
  
  if (is.null(lambda)) stop("lambda can not be empty")
  if (is.character(lambda)) stop("lambda has to be numeric")
  if (lambda<=0) stop("lambda must be positive")  
  
  if (is.null(alpha)) stop("alpha can not be empty")
  if (is.character(alpha)) stop("alpha has to be numeric")
  if (alpha<=0) stop("alpha must be positive for the exponential distribution") 
  
  if (is.null(beta1)) stop("beta1 can not be empty")
  if (is.character(beta1)) stop("beta1 has to be numeric")
  if (is.null(beta2)) stop("beta2 can not be empty")
  if (is.character(beta2)) stop("beta2 has to be numeric")
  if (is.null(beta3)) stop("beta3 can not be empty")
  if (is.character(beta3)) stop("beta3 has to be numeric")
  if (is.null(beta4)) stop("beta4 can not be empty")
  if (is.character(beta4)) stop("beta4 has to be numeric")
  if (is.null(beta5)) stop("beta5 can not be empty")
  if (is.character(beta5)) stop("beta5 has to be numeric")
  if (is.null(gamma)) stop("gamma can not be empty")
  if (is.character(gamma)) stop("gamma has to be numeric")
  
  #Generate covariates following normal distribution with mean 1 and standard deviation 0.5
  V <- rnorm(n, 1, 0.5)
  
  #Generate censoring
  C<-runif(n,min=0,max=max.c)
  
  #Generate indicator for first-stage therapies A1/A2
  X <- rbinom(n, 1, pi.x) # X=0 for A1, 1 for A2
  
  #Generate time to response
  TR <- rexp(n, alpha)
  
  #Generate response status indicator
  ui <- runif(n, 0, 1)
  R <- as.numeric((-log(1-ui)/lambda) >= (TR*exp(beta1*X+gamma*V)))
  
  #Generate indicator for second-stage therapies B1/B2
  Z<-rep(0,n)
  Z[which(R==1)]<-rbinom(length(which(R==1)),1,pi.z)

  #Calculate e1 and e2
  e1 <- exp(beta1*X+gamma*V)
  e2 <- exp(beta1*X+beta2+beta3*X+beta4*Z+beta5*X*Z+gamma*V)
  
  #Generate survival time T
  T<-rep(NA,n)
  T[which(R==0)] <- -log(1-ui[which(R==0)]) / (lambda*e1[which(R==0)])
  T[which(R==1)] <- (-log(1-ui[which(R==1)]) - TR[which(R==1)]*lambda*e1[which(R==1)]) / (lambda*e2[which(R==1)]) + TR[which(R==1)]
  
  #Generate observed time U
  U<-pmin(T,C)
  
  #Generate censoring indicator
  #delta=0 for censoring and delta=1 for event
  delta<-as.numeric(T<C)
  
  #Set time to response to be 0 for non-responders
  TR[which(R==0)] <- 0
  
  #Combine data
  data<-data.frame(X,TR,R,Z,U,delta,V)
  return(data)

}
