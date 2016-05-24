

gen4PMData <- function(NSubj = NULL, abcdParams, D = 1.702,  seed=NULL, theta = NULL, thetaMN = 0, thetaVar = 1){
 ##  Date: January 23, 2016
 ##  Author: Niels Waller
 ##  A program for generating binary item responses for the 1, 2, 3 & 4PM
 ##  
 ##  NSubj       :   Desired number of simulated item response vectors 
 ##  abcdParams  :   a Nitems x 4 matrix of item parameters
 ##  D           :   Deafult = 1.702. Scaling constant to place IRF on normal ogive metric
 ##  seed        :   optional seed for theta generation
 ##  theta       :   User-supplied vector of latent trait scores.  If theta = NULL then
 ##                  gen4PMData with draw NSubj values from a Normal( mean = thetaMN, var = thetaVar) 
 ##                  distribution
  
  
if(!is.null(theta)) NSubj <- length(theta)

## if theta = NULL generate NSubj random normal deviates
  if(is.null(theta)){
      if(is.null(seed)) seed <- sample(1:1000, 1)
      set.seed(seed)
      theta <- sort(rnorm(n=NSubj,mean=thetaMN, sd=sqrt(thetaVar)))
  }    

## Prob of a keyed response for 4PL
P4 <- function(theta,abcd, D){
      a <- abcd[1]
      b <- abcd[2]
      c <- abcd[3]
      d <- abcd[4]
      c + (d-c)/(1+exp(-D*a*(theta-b)))
   }

# create 0/1 observed item response data
# U = probability of a keyed response
# data = 0/1 
NItems <- nrow(abcdParams)
data <- U <- matrix(0,NSubj, NItems)

  for(i in 1:NItems){
      U[,i] <- P4(theta,abcdParams[i,], D)
      data[runif(NSubj) <= U[,i],i] <- 1
  }

 list(data = data, theta = theta, seed = seed)
} ## END gen4PM Data




