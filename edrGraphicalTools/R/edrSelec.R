#Sélection de variable avec SIR
edrSelec <- function(Y, X, H, K, method, pZero=NULL, NZero=NULL, zeta=NULL,
 rho=NULL, baseEst=NULL, btspSamp=NULL, lassoParam=NULL) {

	#Vérification des paramètres d'entrée
	dimVar <- .check.param(Y, X, H, K, "SIR-I")
	n <- dimVar$n
	p <- dimVar$p
	
	if (!(method %in% c("CSS", "RSIR", "SR-SIR"))) {
		stop("Unknown method for selecting components of 'X' with 'SIR-I'.")
	}
	if (all(method=="CSS", any( length(pZero)!=1, length(NZero)!=1))) {
		stop("Invalid imputs 'pZero' and/or 'NZero' for the 'CSS' method.") 
	}
	if (all(method=="CSS", length(zeta)!=1, length(rho)!=1)) {
		stop("Using the 'CSS' method, a value is required for 'zeta' or 'rho'.")
	}
	if (all(method=="CSS", class(baseEst)!="edr", !is.null(baseEst))) {
		stop(paste("Using the 'CSS' method, if an object 'baseEst' is provided,", 
			"it should be an 'edr' object.")) 
	}
	if (all(method=="CSS", K>1)) {
		stop("The 'CSS' method does not work for 'K>1'.")
	}
	if (all(method=="RSIR", baseEst$method!="SIR-I, RSIR", n <= p)) {
		stop(paste("Invalid input 'baseEst' for the RSIR method when n < p.", 
			"Please run the 'edrUnderdet' function."))
	}
	if (all(method=="RSIR", length(btspSamp)!=1)) {
		stop("Invalid input 'btspSamp' for the RSIR method.")
	}
	if (all(method=="SR-SIR", baseEst$method!="SIR-I, SR-SIR", n <= p)) {
		stop(paste("Invalid input 'baseEst' for the SR-SIR method when n < p.", 
			"Please run the 'edrUnderdet' function."))
	}
	if (all(method=="SR-SIR", is.null(lassoParam))) {
		stop("Invalid input 'lassoParam' for the SR-SIR method.")
	}
	
	additionalContent <- NULL
	if (is.null(baseEst)) {
		baseEst <- edrUnderdet(Y, X, H, K, "SIR-QZ")
	}

	### Méthode CSS ###
	if (method=="CSS") {
		if (is.null(baseEst$indices)) {
			baseEst$indices <- X %*% baseEst$matEDR[,1:K]
		}
		#Evaluation des sous-modèles
		matModels <-matrix(0,ncol= pZero,nrow=NZero)
		vectSqCor <- vector("numeric", NZero)
		for (i in 1:NZero){
			matModels[i,] <- sample(1:p, size=pZero,replace=FALSE)
			currEst <- edrUnderdet(Y,X[,matModels[i,]],H,K,"SIR-QZ")
			if (is.null(currEst$indices)) {
				currEst$indices <- X[,matModels[i,]] %*% currEst$matEDR[,1:K]
			}
			vectSqCor[i]<-cor(currEst$indices, baseEst$indices)^2
		}
		#Sélection des meilleurs modèles
		if (is.null(zeta)) {
			matModTop <- matModels[vectSqCor > rho,]
		} else {
			matModTop <- matModels[
				order(vectSqCor,decreasing = TRUE)[1:floor(zeta*NZero)],]
		}	
		scoreVar <- tabulate(matModTop, p)#/tabulate(matIndice.ordonnee)
		additionalContent <- list(matModels=matModels, 
			matModTop=matModTop, vectSqCor=vectSqCor)
	}

	### Méthode RSIR ###
	if (method=="RSIR") {
		if (is.null(baseEst)) {
			baseEst$s <- 0
		}
		#Calcul des échantillons bootstrap
		btspEst <- array(dim=c(p, btspSamp, K))
		for (i in 1:btspSamp) {
			index <- sample(n, n, replace=TRUE)
			M <- sliceMat(Y[index], X[index,], H)
			SigmaTilde <- var(X[index,]) + baseEst$s * diag(p)
			temp <- .edrNorm(M, SigmaTilde,K)
			btspEst[,i,] <- temp$matEDR
		}
		#Réorientation des directions EDR estimées
		baseDir <- (var(X) + baseEst$s * diag(p)) %*% baseEst$matEDR		
		for (k in 1:K) {
    	dd <- (t(baseDir[,k]) %*% btspEst[,,k] > 0) * 2 - 1
    	btspEst[,,k] <- btspEst[,,k] *  as.matrix(rep(1, p)) %*% dd
		}
		#Calcul de 1 - p-value
		scoreVar <- vector("numeric", p)
		for (i in 1:p) {
			betacov <- var(btspEst[i,,])
    	stats <- t(baseEst$matEDR[i,1:K]) %*% solve(betacov) %*% 
				as.matrix(baseEst$matEDR[i,1:K])
    	scoreVar[i] <- pchisq(stats, K)
	}	}

	### Méthode SR-SIR ###
	if (method=="SR-SIR") {
		Sigma <- var(X)
		temp <- sliceMat(Y, X, H, details=TRUE)
		Xh <- t(temp$Xh)
		Ph <- temp$Ph
	
		#Fabrication des éléments de la régression Lasso
		yStar <- Xh %*% sqrt(Ph)
		B <- baseEst$matEDR[,1:K]
		C4 <- B %*% solve(t(B) %*% Sigma %*% Sigma %*% B) %*% t(B) %*% Sigma %*% yStar
		xStar <- NULL
		for (h in 1:H) {
			xStar <- rbind(xStar, Sigma %*% diag(C4[,h]))
		}
		yStar <- as.vector(yStar)

		#Régression Lasso
		aicVec <- vector("numeric", length(lassoParam))
		bicVec <- vector("numeric", length(lassoParam))
		ricVec <- vector("numeric", length(lassoParam))
		matEDR <- list()
		alphaHat <- matrix(0, ncol=length(lassoParam), nrow=p)
		for (i in 1:length(lassoParam)) {
			resL1ce<-l1ce(yStar ~ xStar - 1, sweep.out=NULL, standardize=FALSE, 
				bound=lassoParam[i], absolute.t=TRUE)
   		alphaHat[,i] <- coef(resL1ce)
			matEDR[[i]] <- diag(alphaHat[,i]) %*% baseEst$matEDR[,1:K]
   		rss<-sum((yStar - fitted(resL1ce))^2)


  	 	# Critères de sélection du paramètre Lasso
	   	nVar <-sum(alphaHat[,i] != 0)
  	 	aicVec[i] <- p * H * log(rss / p * H) + 2 * nVar
  	 	bicVec[i] <- p * H * log(rss / p * H) + log(p * H) * nVar
  		ricVec[i] <- (p * H - nVar) * log(rss / (p * H - nVar)) + 
				nVar * (log(p * H) - 1) +	4 / (p * H - nVar - 2)
		}

		chosenPar <- which(ricVec == min(ricVec))[1]
		scoreVar <-	alphaHat[,chosenPar] != 0
		additionalContent <- list(aic=aicVec,bic=bicVec,ric=ricVec,matEDR=matEDR)

	}	


	res <-c(list(scoreVar=scoreVar, K=K, H=H, n=n, method=method, X=X, Y=Y),
		additionalContent)
	class(res) <- "edrSelec"
	res

}

#Fonctions d'affichage pour des objets de type 'edrSelec'
plot.edrSelec <- function(x, nVar=25, ...) {

	if (!inherits(x, "edrSelec")) {
		stop("Use only with 'edrSelec' objects.")
	}
	
	nVar <- min(nVar, length(x$scoreVar))
	if (x$method=="SR-SIR") {
		toShow <- vector("numeric", length(x$scoreVar))
		alphaHat <- NULL
		for (i in 1:length(x$matEDR)) {
			alphaHat <- cbind(alphaHat,diag(x$matEDR[[i]]%*%t(x$matEDR[[i]])) != 0)
		}
		for (i in 1:length(x$scoreVar)) {
			toShow[i] <- min(x$ric) - min(x$ric[alphaHat[i,]], Inf)
		}		
		yLabel <- "Loss in the RIC criterion to select it"
	} else {
		toShow <- x$scoreVar
		if (x$method=="CSS") {
			yLabel <- "Presence in the best models"
		}
		if (x$method=="RSIR") {
			yLabel <- "1 - 'p-value'"
		}
	}
	index <- sort(toShow, decreasing=TRUE, index.return=TRUE)
	par(mar=c(5,5,5,1))
	title <- paste("Variable selection using the", x$method, "method")
	plot(1:nVar, index$x[1:nVar], type="n", xlab="Variable", ylab=yLabel,
		cex.axis=1.5, cex.lab=1.5, cex.main=1.5, main=title, xaxt="n")
	text(1:nVar, index$x[1:nVar], labels=index$ix[1:nVar])

}
summary.edrSelec <- function(object, nVar=5, ...) {

	if (!inherits(object, "edrSelec")) {
		stop("Use only with 'edrSelec' objects.")
	}

	cat(paste("Selection method performed:", object$method),"\n")

	cat(paste("Number of observations:", object$n),"\n")
	cat(paste("Dimension reduction K:", object$K),"\n")
	cat(paste("Number of slices:", paste(object$H, collapse=", ")),"\n")
	cat(" \n")

	if(object$method == "SR-SIR") {
		cat(paste("Selected variables:\n",paste(which(object$scoreVar), 
			collapse=", "),"\n"))
	} else {
		nVar <- min(nVar, length(object$scoreVar))
		cat(paste(nVar,"most important variables:\n"))
		index <- sort(object$scoreVar, decreasing=TRUE, index.return=TRUE)
		toShow <- cbind(index$ix, index$x)
		dimnames(toShow)[[1]] <- 1:length(object$scoreVar)
		dimnames(toShow)[[2]] <- c("Variable", "Score")
		toShow <- toShow[1:nVar,]
		prmatrix(signif(toShow,3))
		cat("\n")
}	}

print.edrSelec <- function(x, ...) {

	if (!inherits(x, "edrSelec")) {
		stop("Use only with 'edrSelec' objects.")
	}
	
	cat(paste("Selection method performed:", x$method),"\n")

	cat(paste("Number of observations:", x$n),"\n")
	cat(paste("Dimension reduction K:", x$K),"\n")
	cat(paste("Number of slices:", paste(x$H, collapse=", ")),"\n")
	cat(" \n")

	if(x$method == "SR-SIR") {
		cat(paste("Selected variables:\n", paste(which(x$scoreVar), 
			collapse=", "),"\n"))
	} else {
		cat("Score of each variable:\n")
		index <- sort(x$scoreVar, decreasing=TRUE, index.return=TRUE)
		toShow <- cbind(index$ix, index$x)
		dimnames(toShow)[[1]] <- 1:length(x$scoreVar)
		dimnames(toShow)[[2]] <- c("Variable", "Score")
		prmatrix(signif(toShow,3))
		cat("\n")
}	}



