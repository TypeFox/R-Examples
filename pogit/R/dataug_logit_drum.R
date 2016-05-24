#### (invisible function)
#### last change: 2015/03/04
#### part I of data_augmentation for binomial logit models: 
#### rewriting the binomial logit model as a latent difference random utility model (dRUM)
#### as in Fussl et al. (2013) based on package "binomlogit"

dataug_binom_dRUM1 <- function(y, N){
  # mixture components (data augmentation)
  # compmix.R taken from the binomlogit package
  
  Hc <- 5 # maximum number of components 
  n  <- length(y)
  ind1 <- as.numeric(y > 0)
  ind2 <- as.numeric(y < N)
  
  mixture <- cbind(matrix(0, nrow = n, ncol = 5), matrix(1, nrow = n, ncol = 5))
  for(s in seq_len(n)){
    mixture[s, 1:length(compmix(N[s])$probs)] <- compmix(N[s])$probs
    mixture[s, 6:(6 + length(compmix(N[s])$probs)-1)] <- compmix(N[s])$var
  }
  probs   <- mixture[, 1:5]
  vars    <- mixture[, 6:10]
  vecvars <- as.vector(t(vars)) 
  
  st <- Hc*(0:(n - 1))
  sd <- sqrt(vars)
  pos   <- numeric(n)
  sr2   <- numeric(n)
  trickmat <- outer(seq_len(Hc), seq_len(Hc), "<=")
  vertfkt  <- matrix(0, nrow = n, ncol = Hc)
  pre      <- probs/sd
  
  # returns list-object containing paramters of the components for normal mixture  
  comp.aug <- list(Hc = Hc, pre = pre, vars = vars, vecvars = vecvars,
                   trickmat = trickmat, st = st, ind1 = ind1, ind2 = ind2)
  return(comp.aug)
}


#### (invisible function)
#### last change: 2015/03/04
#### part II of data_augmentation for binomial logit models: 
#### rewriting the binomial logit model as a latent difference random utility model (dRUM)
#### as in Fussl et al. (2013) based on package "binomlogit"

dataug_binom_dRUM2 <- function(y, N, mu, comp.aug){
	n <- length(y)
	lambda <- exp(mu)
		
	u <- rgamma(n, shape = N, rate = 1 + lambda)
	v <- rgamma(n, shape = y, rate = 1)
	w <- rgamma(n, shape = N - y, rate = lambda)
	ystar <- -log((u + comp.aug$ind2*w)/(u + comp.aug$ind1*v))
	
	# non-standardized probabilities 	
	Pr <- comp.aug$pre * exp(rep.int((ystar - log(lambda))^2, comp.aug$Hc)/-2*comp.aug$vars)
	# non-standardized cdfs (=> inversion method)
	vertfkt <- Pr %*% comp.aug$trickmat
	pos     <- rowSums(vertfkt[, comp.aug$Hc]*runif(n) > vertfkt) + 1
	sr2     <- comp.aug$vecvars[comp.aug$st + pos]
	invSig  <- 1/sqrt(sr2)
	return(list(ystar = ystar, invSig = invSig)) 
}




