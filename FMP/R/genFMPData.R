

genFMPData<-function(NSubj, bParams, theta = NULL, thetaMN = 0, thetaVar = 1, seed){

###################################################################
# GENERATE Data FOR testing FMP routines
# SAME DATA USED IN BROWNE AND WALLER N=2000 CYN MMPI-A EXAMPLE
  
  if(is.data.frame(bParams)) bParams <- as.matrix(bParams)

  if(ncol(bParams)!=9) stop("\n\nbParams should have 9 columns")

  # Prob of a keyed response 2PL ------------------------------##
  P <- function(m){
     1/(1+exp(-m))
  }

  NItems = nrow(bParams)
  set.seed(seed)

  # if theta not supplied by the user, generate NSubj random 
  # normal deviates from a population with mean = thetaMN and
  # variance = thetaVar
  if(is.null(theta)){
       x <- sort(rnorm(n=NSubj,mean=thetaMN, sd=thetaVar))
  }
  else x <- theta

  xpoly <- matrix(cbind(1,x, x^2, x^3, x^4, x^5, x^6, x^7),nrow=NSubj, 8)


  # create 0/1 observed item response data
  # U = probability of a keyed response
  # data = 0/1 
  data <- U <- matrix(0,NSubj, NItems)
  for(item in 1:NItems){
     U[,item] <- P(xpoly %*% bParams[item,1:8])
     data[runif(NSubj) <=U[,item],item] <- 1
  }

  list(theta = x, data = data, seed = seed)

} ## END genFMPData

