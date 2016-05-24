drs <- function(formula, data, type = c("d", "rs"), C = NULL,	pooled = FALSE) {
	
  call <- match.call()
  type <- match.arg(type)
	dat <- data
	mf <- model.frame(formula = formula, data = dat)
	resp <- model.response(mf)
	m <- resp[, 1]
	n <- rowSums(resp)
  dat <- data.frame(n = n, m = m, group = factor(mf[, 2]))
  rho <- NULL
	
	#--- donner
	if(type == "d") {
    # computation of rho and correction factors C
  	group <- sort(unique(dat$group))
  	lev <- levels(group)
  	N <- ntot <- mtot <- mu <- nA <- 0
  	dat$mu <- rep(0, nrow(dat))
  	for(i in seq(length(group))){
    	tmp <- dat[dat$group == lev[i], ]
    	n <- tmp$n ; m <- tmp$m
    	N[i] <- nrow(tmp)
  		ntot[i] <- sum(n)
  		mtot[i] <- sum(m)
  		mu[i] <- mtot[i] / ntot[i]
    	dat$mu[dat$group == lev[i]] <- mu[i]
    	nA[i] <- sum(dat$n[dat$group == lev[i]]^2) / ntot[i]
  	}
  	df.SSC <- sum(N) - length(group)
  	df.SSE <- sum(ntot) - sum(N)
  	MSC <- sum(dat$n * (dat$m / dat$n - dat$mu)^2) / df.SSC
  	MSE <- sum(dat$n * (dat$m / dat$n) * (1 - dat$m / dat$n)) / df.SSE
  	K <- (sum(ntot) - sum(nA)) / (sum(N) - length(group))
  	rho <- (MSC - MSE) / (MSC + (K - 1) * MSE)
  	if(is.null(C))
    	C <- 1 + (nA - 1) * rho
  
		# results
  	tab <- data.frame(group = group, N = N, n = ntot, m = mtot, mu = mu, C = C)   
  	names(tab)[1] <- as.character(attr(terms(formula), "variables"))[3]
  	mu <- sum(tab$m) / sum(tab$n)
  	X2 <- sum((tab$m - tab$n * mu)^2 / (tab$C * tab$n * mu * (1 - mu)))
	}
	
	#--- rs
	if(type == "rs") {
		# computation of mu, vratio and vbin per group
	  group <- sort(unique(dat$group))
  	lev <- levels(group)
  	N <- ntot <- mtot <- mu <- vbin <- vratio <- 0
  	for(i in 1:length(group)){
    	tmp <- dat[dat$group == lev[i], ]
    	n <- tmp$n ; m <- tmp$m
    	N[i] <- nrow(tmp) ; ntot[i] <- sum(n) ; mtot[i] <- sum(m)
    	mu[i] <- mtot[i] / ntot[i] ;
    	vbin[i] <- mu[i] * (1 - mu[i]) / ntot[i]
    	vratio[i] <- N[i] * (N[i] - 1)^(-1) * ntot[i]^(-2) * sum((m - n * mu[i])^2)
  	}

		# computation of the design effect C
  	if(is.null(C)){
    	C <- vratio / vbin
    	if(pooled){
      	df <- length(group) - 1
      	mumean <- sum(mtot) / sum(ntot)
      	A <- (1 / df) * (1 / (mumean * (1 - mumean)))
      	B <- (1 - ntot / sum(ntot)) * mu * (1 - mu) * C
      	d <- A * sum(B) ; C <- rep(d, length(C))
    	}
  	}

		# results
  	tab <- data.frame(group = group, N = N, n = ntot, m = mtot, mu = mu, 
  		vbin = vbin, vratio = vratio, deff = C)
  	names(tab)[1] <- as.character(formula[3])
  	nadj <- tab$n / tab$deff
		madj <- tab$m / tab$deff
  	muadj <- sum(madj) / sum(nadj)
  	X2 <- sum((madj - nadj * muadj)^2 / (nadj * muadj * (1 - muadj)))
		
	}

  # outputs
  structure(
  	list(call = call, type = type, tab = tab, rho = rho, X2 = X2, dat = dat),
  	class = "drs"
  	)
  
}
