payoff <- function(payType, pcFlag, strike, S) {
	
  # check inputs
  if (nargs() != 4) {
    stop("incorrect number of arguments")
  }
  if(!is.numeric(pcFlag) | !is.numeric(payType) | 
       !is.numeric(strike)) {
    stop("pcFlag, payType, and strike must be numeric")
  }
  if (length(payType) != 1) {
    stop("length of payType must be 1")
  }
  if(is.list(S) == FALSE) {
    stop(paste("S must be a list containing the vectors of the spatial",
        "grid points of the underlying assets"))
  } else {
    # determine number of underlying assets
    nAsset <- length(S)
  }
  if(length(pcFlag) != nAsset | length(strike) != nAsset) {
    stop(paste("length of pcFlag and strike must equal the",
               "number of underlying assets"))
  }

  # find dimension of underlying spatial nodes
  dimS <- c()
  for (i in 1:nAsset) {
    dimS <- c(dimS, length(S[[i]]))
  }
  
	# calculate option payoff
  
  # digital
	if (payType == 0) {

    moneyness <- c()
  	for (i in 1:nAsset) {
  		if (pcFlag[i] == 0) {
  			moneyness <- c(moneyness, list(S[[i]] >= strike[i]))
  		} else if (pcFlag[i] == 1) {
  			moneyness <- c(moneyness, list(S[[i]] <= strike[i]))
  		} else {
  			stop("invalid pcFlag. Must be 0 or 1")
  		}
  	}
  	outerMoneyness <- array(apply(expand.grid(moneyness), 1, sum), 
                            dim=dimS)
  	payoff <- (outerMoneyness == nAsset) * 1

	# best-of
	} else if (payType == 1) {
    
	  moneyness <- c()
	  for (i in 1:nAsset) {
	    if (pcFlag[i] == 0) {
	      moneyness <- c(moneyness, list(pmax(S[[i]] - strike[i], 0)))
	    } else if (pcFlag[i] == 1) {
	      moneyness <- c(moneyness, list(pmax(strike[i] - S[[i]], 0)))
	    } else {
	      stop("invalid pcFlag. Must be 0 or 1")
	    }
	  }
	  payoff <- array(apply(expand.grid(moneyness), 1, max), dim=dimS)

	# worst-of
	} else if (payType == 2) {
    
	  moneyness <- c()
	  for (i in 1:nAsset) {
	    if (pcFlag[i] == 0) {
	      moneyness <- c(moneyness, list(pmax(S[[i]] - strike[i], 0)))
	    } else if (pcFlag[i] == 1) {
	      moneyness <- c(moneyness, list(pmax(strike[i] - S[[i]], 0)))
	    } else {
	      stop("invalid pcFlag. Must be 0 or 1")
	    }
	  }
	  payoff <- array(apply(expand.grid(moneyness), 1, min), dim=dimS)

	# payType error catch
	} else {
		stop("invalid payType. Must be 0, 1, or 2")
	}

	# return payoff array
	payoff

}