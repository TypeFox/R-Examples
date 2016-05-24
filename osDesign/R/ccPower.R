ccPower <-
function(B=1000,
										betaTruth,
										X,
										N,
  	 	              expandX="all",
	   	              etaTerms=NULL,									 
										nCC,
										r=1,
										alpha=0.05,
                    digits=1,
                    betaNames=NULL,
                    monitor=NULL)
{
	##
	problem <- coreChecks(betaTruth=betaTruth, X=X, N=N, etaTerms=etaTerms, expandX=expandX, betaNames=betaNames)
	if(problem != "")
		stop(problem)

	##
	problem <- ccChecks(nCC=nCC)
	if(problem != "")
		stop(problem)

  ##
  if(is.null(colnames(X)))
  	colnames(X) <- c("Int", paste("V", 1:(ncol(X) - 1), sep = ""))
	if(is.null(monitor))
		monitor <- B + 1
  if(length(r) > 1)
		cat("\nWarning: taking the first value in 'r'\n")

	## Restrict design matrix to columns indicated in 'etaTerms'
	##
	if(!is.null(etaTerms))
		X <- X[, is.element(colnames(X), etaTerms)]

	## Restrict design matrix to columns indicated in 'etaTerms'
	##
	if(!is.null(etaTerms))
		X <- X[, is.element(colnames(X), etaTerms)]

	## Expanded design matrix
  ##
  Xexp <- X
  if(expandX[1] != "none")
		Xexp <- as.data.frame(expandCatX(X, expandX=expandX))
	if(is.null(betaNames))
		betaNames <- names(Xexp)

	## Outcome probabilities
  ##
  piY <- expit(as.vector(as.matrix(Xexp) %*% betaTruth))
  
  ## Formulae for glm() calls
  ##
	formCD <- as.formula(paste("cbind(N1, N0) ~", paste(colnames(Xexp)[-1], collapse=" + ", sep="")))
	formCC <- as.formula(paste("cbind(n1, n0) ~", paste(colnames(Xexp)[-1], collapse=" + ", sep="")))
	
	## Run simulation
	##
  nDesigns <- 1 + length(nCC)
  p        <- length(betaTruth)
	waldTest <- array(NA, dim=c(B, nDesigns, p))
  cat(paste(nDesigns, "designs will be simulated\n"))
	##
	for(b in 1:B)
	{
		##
    if((b %% monitor) == 0)
      cat("Repetition", b, "of", B, "complete\n")
    
		##
		Xexp$N1 <- rbinom(nrow(Xexp), N, piY)
		Xexp$N0 <- N - Xexp$N1
		fitCD   <- summary(glm(formCD, data=Xexp, family=binomial()))$coef
		waldTest[b,1,] <- (fitCD[,4] < alpha)
		
		##
		for(i in 1:length(nCC))
		{
			n0 <- round(nCC[i] * r[1] / (r[1]+1))
			n1 <- nCC[i] - n0
			Xexp$n0 <- rmvhyper(Xexp$N0, n0)
			Xexp$n1 <- rmvhyper(Xexp$N1, n1)
	    XexpCC  <- Xexp[(Xexp$n0 > 0 | Xexp$n1 > 0),]
			##
			fitCC   <- try(glm(formCC, data=XexpCC, family=binomial()), silent=TRUE)
			if(class(fitCC)[1] == "glm")
			{
				fitCC <- summary(fitCC)$coef
				if(nrow(fitCC) == p)
					waldTest[b,(1+i),] <- c(FALSE, (fitCC[-1,4] < alpha))
			}
		}
	}
	
	## Error checks for output and evaluate power
	##
	keep    <- apply(is.na(waldTest), c(1,2), sum) == 0
	results <- evalPower(waldTest, keep, resNames=list(c("CD", paste("CC n =", nCC, "  ")), betaNames))
	failed  <- B - matrix(apply(keep, 2, sum), ncol=1, dimnames=list(rownames(results), ""))
	
	
	## Return object of class 'ccPower'
  ##
  value           <- NULL
  value$B         <- B
  value$betaTruth <- betaTruth
  value$X         <- X
  value$N         <- N
  value$nCC       <- nCC
  value$r         <- r
  value$alpha     <- alpha
	value$digits 	  <- digits
  ##
  value$failed    <- failed
  value$betaPower <- results
  ##
  class(value) <- "ccPower"
  return(value)
}
