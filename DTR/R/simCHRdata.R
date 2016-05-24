###################################################
### Reference:
### Tang X, Wahed AS: Cumulative hazard ratio estimation for treatment regimes in
### sequentially randomized clinical trials. Statistics in Biosciences, [Epub ahead of print]
###################################################

###################################################
### code chunk number 1: chunklibraries
###################################################

#Libraries required
require(survival)

###################################################
### code chunk number 2: chunksimulatedata
###################################################

simCHRdata <- function(n, # Number of subjects in the clinical trial
                         max.c, # Censoring time C is generated from Uniform(min=max.c/2,max=max.c)
                         pi.x, # First-stage indicator X is generated from Bernoulli(pi.x); X=0 for assignment to A1 and X=1 for assignment to A2
                         pi.r, # Remission/Consent indicator R is Bernoulli(pi.r)
                         pi.z, # Second-stage indicator Z is generated from Bernoulli(pi.z); Z=0 for assignment to B1 and Z=1 for assignment to B2
                         gamma10, # When X=0, R=0, survival time is draw from weibull with alpha10 and gamma10
                         gamma11, # When X=0, R=1, Z=0, survival time is drawn from weibull with alpha11 and gamma11
                         gamma12, # When X=0, R=1, Z=1, survival time is drawn from weibull with alpha12 and gamma12
                         gamma20, # When X=1, R=0, survival time is draw from weibull with alpha20 and gamma20
                         gamma21, # When X=1, R=1, Z=0, survival time is drawn from weibull with alpha21 and gamma21
                         gamma22, # When X=1, R=1, Z=1, survival time is drawn from weibull with alpha22 and gamma22
                         alpha10, # When X=0, R=0, survival time is draw from weibull with alpha10 and gamma10
                         alpha11, # When X=0, R=1, Z=0, survival time is drawn from weibull with alpha11 and gamma11
                         alpha12, # When X=0, R=1, Z=1, survival time is drawn from weibull with alpha12 and gamma12
                         alpha20, # When X=1, R=0, survival time is draw from weibull with alpha20 and gamma20
                         alpha21, # When X=1, R=1, Z=0, survival time is drawn from weibull with alpha21 and gamma21
                         alpha22, # When X=1, R=1, Z=1, survival time is drawn from weibull with alpha22 and gamma22
                         beta # Coefficiant vector for two covariates, for example, beta=c(0.5,0.5)
) {

  #Functions needed to check errors:
  #Function to test if a number is an integer or not.
  is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  #Function to check if a number is a valid probability value or not
  is.probability<-function(x) if(x>=0 && x<=1) TRUE else FALSE
  
  #Check input data
  if (is.null(n)) stop("n can not be empty")
  if (is.character(n)) stop("n has to be numeric")
  if (n<=0 || !is.wholenumber(n)) stop("n must be a positive integer")
  
  if (is.null(max.c)) stop("max.c can not be empty")
  if (is.character(max.c)) stop("max.c has to be numeric")
  if (max.c<=0) stop("max.c must be a positive value")  
  
  if (is.null(pi.x)) stop("pi.x can not be empty")
  if (is.character(pi.x)) stop("pi.x has to be numeric")
  if (!is.probability(pi.x)) stop("pi.x must be a valid probability value") 
  
  if (is.null(pi.r)) stop("pi.r can not be empty")
  if (is.character(pi.r)) stop("pi.r has to be numeric")
  if (!is.probability(pi.r)) stop("pi.r must be a valid probability value")   
  
  if (is.null(pi.z)) stop("pi.z can not be empty")
  if (is.character(pi.z)) stop("pi.z has to be numeric")
  if (!is.probability(pi.z)) stop("pi.z must be a valid probability value")
  
  if (is.null(gamma10)) stop("gamma10 can not be empty")
  if (is.character(gamma10)) stop("gamma10 has to be numeric")
  if (is.null(gamma11)) stop("gamma11 can not be empty")
  if (is.character(gamma11)) stop("gamma11 has to be numeric")
  if (is.null(gamma12)) stop("gamma12 can not be empty")
  if (is.character(gamma12)) stop("gamma12 has to be numeric")
  if (is.null(gamma20)) stop("gamma20 can not be empty")
  if (is.character(gamma20)) stop("gamma20 has to be numeric")
  if (is.null(gamma21)) stop("gamma21 can not be empty")
  if (is.character(gamma21)) stop("gamma21 has to be numeric")
  if (is.null(gamma22)) stop("gamma22 can not be empty")
  if (is.character(gamma22)) stop("gamma22 has to be numeric")
  
  if (is.null(alpha10)) stop("alpha10 can not be empty")
  if (is.character(alpha10)) stop("alpha10 has to be numeric")
  if (is.null(alpha11)) stop("alpha11 can not be empty")
  if (is.character(alpha11)) stop("alpha11 has to be numeric")
  if (is.null(alpha12)) stop("alpha12 can not be empty")
  if (is.character(alpha12)) stop("alpha12 has to be numeric")
  if (is.null(alpha20)) stop("alpha20 can not be empty")
  if (is.character(alpha20)) stop("alpha20 has to be numeric")
  if (is.null(alpha21)) stop("alpha21 can not be empty")
  if (is.character(alpha21)) stop("alpha21 has to be numeric")
  if (is.null(alpha22)) stop("alpha22 can not be empty")
  if (is.character(alpha22)) stop("alpha22 has to be numeric")  
  #Check beta
  if (is.null(beta)) stop("beta can not be empty")
  for(i in 1:length(beta)) { 
    if (is.character(beta[i])) stop("beta has to be numeric")
  }
  
  #Generate two covariates following the bernoulli(0.5)
  V1 <- rbinom(n, 1, 0.5)
  V2 <- rbinom(n, 1, 0.5)
  
  #Generate censoring time
  C <- runif(n, min=max.c/2, max=max.c)

  #Generate indicator for first-stage therapies A1/A2
  X <- rbinom(n, 1, pi.x) # X=0 for A1, 1 for A2
  
  #Generate response indicator R=0 nonresponders R=1 responders
  R <- rbinom(n, 1, pi.r)
  
  #Generate indicator for second-stage therapies B1/B2
  Z<-rep(0,n)
  Z[which(R==1)]<-rbinom(length(which(R==1)),1,pi.z) # Z=0 for B1, 1 for B2
  
  #Calculate exp(beta*V)
  e <- exp(beta[1]*V1 + beta[2]*V2)
  
  #Generate survival time
  ui <- runif(n, min=0, max=1)
  T <- rep(NA, n)
  
  # Generate survival time for A1NR
  T[which(X==0 & R==0)] <- (-log(ui[which(X==0 & R==0)]) / (alpha10 * e[which(X==0 & R==0)]))^(1/gamma10)

  # Generate survival time for A1RB1
  T[which(X==0 & R==1 & Z==0)] <- (-log(ui[which(X==0 & R==1 & Z==0)]) / (alpha11 * e[which(X==0 & R==1 & Z==0)]))^(1/gamma11)

  # Generate survival time for A1RB2
  T[which(X==0 & R==1 & Z==1)] <- (-log(ui[which(X==0 & R==1 & Z==1)]) / (alpha12 * e[which(X==0 & R==1 & Z==1)]))^(1/gamma12)

  # Generate survival time for A2NR
  T[which(X==1 & R==0)] <- (-log(ui[which(X==1 & R==0)]) / (alpha20 * e[which(X==1 & R==0)]))^(1/gamma20)

  # Generate survival time for A2RB1
  T[which(X==1 & R==1 & Z==0)] <- (-log(ui[which(X==1 & R==1 & Z==0)]) / (alpha21 * e[which(X==1 & R==1 & Z==0)]))^(1/gamma21)

  # Generate survival time for A2RB2
  T[which(X==1 & R==1 & Z==1)] <- (-log(ui[which(X==1 & R==1 & Z==1)]) / (alpha22 * e[which(X==1 & R==1 & Z==1)]))^(1/gamma22)

  # Generate observed time U
  U <- pmin(T, C)

  #Generate censoring indicator
  #delta=0 for censoring and delta=1 for event
  delta <- as.numeric(T <= C)

  # Combine data
  data <- data.frame(X, R, Z, U, delta, V1, V2)
  return(data)

}

