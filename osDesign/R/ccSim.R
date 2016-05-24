ccSim <-
function(B=1000,
									betaTruth,
									X,
									N,
   	              expandX="all",
   	              etaTerms=NULL,									 
									nCC,
									r,
									refDesign=1,
									alpha=0.05,
									threshold=c(-Inf, Inf),
                  digits=1,
                  betaNames=NULL,
                  monitor=NULL,
                  returnRaw=FALSE)
{
	##
	problem <- coreChecks(betaTruth=betaTruth, X=X, N=N, etaTerms=etaTerms, expandX=expandX, betaNames=betaNames)
	if(problem != "")
		stop(problem)

	##
	problem <- ccChecks(nCC=nCC, threshold=threshold)
	if(problem != "")
		stop(problem)

  ##
  if(is.null(colnames(X)))
  	colnames(X) <- c("Int", paste("V", 1:(ncol(X) - 1), sep = ""))
  if(length(nCC) > 1)
		cat("\nWarning: taking the first value in 'nCC'")
	if(is.null(monitor))
		monitor <- B + 1
	
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
	rStar    <- unique(c(refDesign, r))
	nDesigns <- length(rStar) + 1
	p        <- length(betaTruth)	
	betaHat  <- array(NA, dim=c(B, nDesigns, p))
	seHat    <- array(NA, dim=c(B, nDesigns, p))
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
		betaHat[b,1,]  <- fitCD[,1]
		seHat[b,1,]    <- fitCD[,2]
		waldTest[b,1,] <- (fitCD[,4] < alpha)
		
		##
		for(i in 1:length(rStar))
		{
			n0 <- round(nCC[1] * rStar[i] / (rStar[i]+1))
			n1 <- nCC[1] - n0
			Xexp$n0 <- rmvhyper(Xexp$N0, n0)
			Xexp$n1 <- rmvhyper(Xexp$N1, n1)
	    XexpCC  <- Xexp[(Xexp$n0 > 0 | Xexp$n1 > 0),]
			##
			fitCC   <- try(glm(formCC, data=XexpCC, family=binomial()), silent=TRUE)
			if(class(fitCC)[1] == "glm")
			{
				fitCC <- summary(fitCC)$coef
				if(nrow(fitCC) == p)
				{
					betaHat[b,(1+i),-1]  <- fitCC[-1,1]
					seHat[b,(1+i),-1]    <- fitCC[-1,2]
					waldTest[b,(1+i),-1] <- (fitCC[-1,4] < alpha)
				}
			}
		}
	}
	
	## Error checks for output and evaluate operating characteristics
	##
	keep    <- keepOC(betaTruth, betaHat, seHat, threshold)
	results <- evalOC(betaTruth, betaHat, seHat, waldTest, keep, alpha=alpha, resNames=list(c("CD", paste("CC r =", rStar, "")), betaNames))
	failed  <- B - matrix(apply(keep, 2, sum), ncol=1, dimnames=list(rownames(results$betaPower), ""))
	
	## Return object of class 'ccSim'
  ##
  value           <- NULL
  value$B         <- B
  value$betaTruth <- betaTruth
  value$X         <- X
  value$N         <- N
  value$nCC       <- nCC
  value$r         <- r
  value$refDesign <- refDesign
  value$alpha     <- alpha
  value$threshold <- threshold
	value$digits 	  <- digits
  ##
  value$failed    <- failed
  value$results   <- results
  ##
  if(returnRaw == TRUE)
  {
  	value$betaHat <- betaHat
		value$seHat   <- seHat
  }
  class(value) <- "ccSim"
  return(value)
}
