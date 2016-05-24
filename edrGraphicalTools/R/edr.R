edr <-
function(Y,X,H,K,method,submethod="SIR-QZ", ...){

	tol <- 1e-10

	#Vérification des paramètres d'entrée
	dimVar <- .check.param(Y, X, H, K, method)
	n <- dimVar$n
	p <- dimVar$p
	
	if(n<=p) {
		if (method=="SIR-I") {
			warning(paste("The sample size 'n' is not greater than the number of",
				"variables 'p'. We will use a specific method dedicated to this case."))
			res <- edrUnderdet(Y, X, H, K, submethod, ...)
		} else {
			stop(paste("The sample size 'n' is not greater than the number of",
				"variables 'p'. This case is not handled when 'method' is not equal to", 
				"'SIR-I'."))
		}
	} else {

		nbind<-rep(n%/%H,H) # number of subject in each slices
		if((n%/%H)*H!=n){
			ad<-rep(1,H)
			ad[sample(H,H-(n-(n%/%H)*H))]<-0
			nbind<-nbind+ad
		}

		res<-svd(((n-1)/n)*var(X))
		sigmainvsqrt<-res$u%*%diag(1/sqrt(res$d))%*%t(res$v)
		
		matXY<-cbind(Y,X) 
		Xord<-matXY[sort.list(matXY[,1]),-1] #sort X by Y

		M1 <- sliceMat(Y, X, H, rdSup=TRUE)
		M1 <- sigmainvsqrt%*%M1%*%sigmainvsqrt
		M2<-matrix(0,ncol=p,nrow=p)
		Msave<-matrix(0,ncol=p,nrow=p)


		Vbar<-matrix(0,ncol=p,nrow=p)

		for(h in 1:H){
			Vh<-((nbind[h]-1)/nbind[h])*var(Xord[rep(1:H,nbind)==h,])
			Vbar<-Vbar+(nbind[h]/n)*Vh
		}

		for(h in 1:H){
			M2<-M2+(nbind[h]/n)*(sigmainvsqrt%*%(Vh-Vbar)%*%sigmainvsqrt)%*%
				t(sigmainvsqrt%*%(Vh-Vbar)%*%sigmainvsqrt)
			Msave<-Msave+(nbind[h]/n)*(diag(p)-sigmainvsqrt%*%Vh%*%sigmainvsqrt)%*%
				t(diag(p)-sigmainvsqrt%*%Vh%*%sigmainvsqrt)
		}
	
		M<-switch(method, 
			"SIR-I"= M1,
			"SIR-II"= M2,
			"SAVE"= Msave)

		rSIR<-eigen(M, symmetric=TRUE)
		matEDR <- as.matrix(sigmainvsqrt %*% rSIR$vectors)
		if (any( abs(Im(matEDR)) > tol )) {
			warning("Some large imaginary parts were deleted in the complex matrix 'matEDR'.")
		}
		matEDR <- Re(matEDR)

		res <- list(matEDR = matEDR, eigvalEDR = rSIR$values, K=K, H=H, n=n, method=method,
								X=X, Y=Y)
		class(res) <-"edr"  
		res
}	}

#Calcul la matrice de tranchange d'une méthode SIR-I
sliceMat <- function(Y, X, H, details=FALSE, rdSup=FALSE) {

	p <- ncol(X)
	n <- nrow(X)
	index <- sort(Y, index.return=TRUE)$ix
	if (rdSup) {
		nbind<-rep(n%/%H,H) 
		if((n%/%H)*H!=n){
			ad<-rep(1,H)
			ad[sample(H,H-(n-(n%/%H)*H))]<-0
			nbind<-nbind+ad
		}
		group <- rep(0:(H-1), nbind)
	} else {		
		group <- sort(1:n %% H)
	}
	  
	Sh <- matrix(0, ncol=H,nrow=n)
	Ph <- matrix(0, ncol=H, nrow=H)
	for (i in 0:(H-1)) {
	    nh <- sum(group==i)
	    Sh[index[group==i],i+1] <- 1/nh
	    Ph[i+1,i+1] <- nh/n
	}
	Sh <- Sh - 1/n # Sh <- "Ibar_n" %*% S_h
	resultat <- t(X) %*% Sh %*% Ph %*% t(Sh) %*% X

	if (details) {
		list(M=resultat, Xh= t(Sh) %*% X, Ph=Ph)
	} else {
		resultat
}	}


#Estime les directions EDR ou les indices d'un modèle de régression lorsque n < p
#en utilisant SIR-I
edrUnderdet <- function(Y, X, H, K, method, initEDR=NULL, maxIter=NULL,
	regulParam=NULL, sMin=1e-16, sChg=10, btspSamp=NULL) {



	tol <- 1e-6

	#Vérification des paramètres d'entrée
	dimVar <- .check.param(Y, X, H, K, "SIR-I")
	n <- dimVar$n
	p <- dimVar$p
	#if (!(method %in% c("SIR-QZ", "SIR-MP", "RSIR", "SR-SIR"))) {
	if (!(method %in% c("SIR-QZ", "RSIR", "SR-SIR"))) {
		stop("Unknown method for 'SIR-I' when n < p.")
	}
	if (all(method!="SIR-QZ", length(H)>1)) {
		warning("Only the first value of the 'H' vector will be considered.")	
		H <- H[1]
	}
	if (all(method=="SR-SIR", is.null(initEDR))) {
		initEDR <- as.matrix(diag(p)[,1:K])
	}
	if (all(method=="SR-SIR", any(length(maxIter)!=1, dim(initEDR)!=c(p,K), 
		is.null(regulParam) ))) {
		stop(paste("Invalid inputs 'initEDR', 'maxIter' and/or 'regulParam'",
			"for the 'SR-SIR' method."))
	} 
	if (all(method=="RSIR", any(length(btspSamp)!=1, is.null(regulParam) ))) {
		stop(paste("Invalid inputs 'btspSamp' and/or 'regulParam'",
			"for the 'RSIR' method."))
	} 
	if (all(method=="SIR-QZ", any(!is.numeric(sMin), !is.numeric(sChg), sMin <= 0, 
		sChg <= 0))) {
		stop("Invalid inputs 'sMin' and/or 'sChg' for the 'SR-SIR' method.")
	}
	if (all( any(method=="RSIR",method=="SR-SIR"), any(regulParam < 0))) {
		stop("Every element in 'regulParam' must be positive.")
	}

	#Exécution des différentes méthodes
	Sigma <- var(X)
	additionalContent <- NULL
	if (method == "SIR-QZ") {
		sMax <- 1e10

		sVect <- NULL
		eigvalEDR <- NULL
		indices <- NULL
		for (i in 1:length(H)) {
			M <- sliceMat(Y, X, H[i])
			temp <- .Rdggev(M, Sigma)
			s <- sMin
			cache <- is.finite(temp$lambda)
			index <- sort(temp$lambda[cache], index.return=TRUE, 
				decreasing=TRUE)$ix			
			while ( (sum((abs(temp$beta) < tol) & (abs(temp$alphar) < tol)) > 0)
          | (length(temp$alphar) - sum(temp$alphar < tol) < K) ) {
				if (s >= sMax) {
					stop("The 'SIR-QZ' methods fails to find a regularization parameter.")
				}
		  	temp <- .Rdggev(M, Sigma + s * diag(p))
  			s <- s * sChg
   			cache <- is.finite(temp$lambda)
		  	index <- sort(temp$lambda[cache], index.return=TRUE, 
					decreasing=TRUE)$ix
			}
			sVect <- c(sVect, s)
			#Quelques valeurs de 'temp$lambda' peuvent éventuellement être infinies, 
			#elles sont retirées avant de renvoyer ces valeurs propres
			if (is.null(eigvalEDR)) {
				eigvalEDR <- as.matrix(temp$lambda[cache][index])
			}	else {
				diffTaille <- dim(eigvalEDR)[1] - length(temp$lambda[cache][index])
				if (diffTaille > 0) {
					eigvalEDR	<- cbind(eigvalEDR, c(temp$lambda[cache][index], 
						rep(NA,diffTaille)))
				} else {
					eigvalEDR	<- cbind(rbind(eigvalEDR, 
						matrix(NA, nrow=-diffTaille, ncol=dim(eigvalEDR)[2])), 
						temp$lambda[cache][index])
			}	}
			indices <- cbind(indices, X %*% temp$vr[,cache][,index[1:K]])
		}
		
		resACP <- princomp(indices)
		indices <- resACP$scores[,1:K]

		res <- list(indices=as.matrix(indices), eigvalEDR=eigvalEDR, K=K, H=H, n=n,
			method=paste("SIR-I,",method), X=X, Y=Y, s=sVect)

	} else {
		#Non fonctionnel :
		if (method == "SIR-MP") {
			M <- sliceMat(Y, X, H)
			M.plus <- ginv(M)
			lambda.stop <- min(H-1,p)
			#A gérer : H < K.

			Sigma.svd <- svd(Sigma)
			Sigma.half <- Sigma.svd$u %*% diag(sqrt(Sigma.svd$d)) %*% Sigma.svd$v

			decElPro <- eigen(Sigma.half %*% M.plus %*% Sigma.half)
			eta <- decElPro$vectors[,lambda.stop + 1 - 1:K]
			matEDR <- M.plus %*% Sigma.half %*% eta
			matEDR <- Re(matEDR %*% diag(1/sqrt(diag(t(matEDR) %*% matEDR))))
			eigvalEDR <- Re(decElPro$values)

		}	
		if (method == "RSIR") {
			regulParam <- sort(regulParam)
			if (regulParam[1]!=0) {
				regulParam <- c(0, regulParam)
			}
			M <- sliceMat(Y, X, H)

			#Directions EDR initiales, sans régularisation
			temp <- .edrNorm(M, Sigma, K)
			matEDRBase <- temp$matEDR
			eigValBase <- temp$eigVal
	
		
			#Détermination du parmètre de lissage par estimation bootstrap du MSE
			#-> Calcul des estimations des directions EDR, des variances et des moyennes
			vars <- matrix(0, nrow=length(regulParam), ncol=K)
			matEDRMean <- array(dim=c(p,length(regulParam), K))
			for (i in 1:length(regulParam)) {
				matEDRBtsp <- array(dim=c(p,btspSamp,K))
				for (j in 1:btspSamp) {
					index <- sample(n,n,replace=TRUE)
					M2 <- sliceMat(Y[index], X[index,], H)
					Sigma2 <- var(X[index,])
      		temp <- .edrNorm(M2, Sigma2 + regulParam[i] * diag(p), K)
					#print(dim(temp$matEDR))
					matEDRBtsp[,j,] <- temp$matEDR
				}
				matEDRBaseReg <- (Sigma + regulParam[i] * diag(p)) %*% matEDRBase
				for (k in 1:K) {
					 D <- (t(matEDRBaseReg[, k]) %*% matEDRBtsp[, , k] > 0) * 2 - 1
      		 matEDRBtsp[,,k] <- matEDRBtsp[,,k] *  as.matrix(rep(1, p)) %*% D
      		 vars[i,k] <- sum(diag(var(t(matEDRBtsp[,,k]))))
			     matEDRMean[,i,k] <- matEDRBtsp[,,k] %*% rep(1/btspSamp,btspSamp)				
			} }			

			#-> Calcul de l'estimation du MSE
			loss <- vector("numeric", length(regulParam))
			bias <- array(dim=c(p,length(regulParam), K))
			for (k in 1:K) {
    		bias[,,k] <- matEDRMean[,,k] - matEDRMean[,1,k] %*% 
					t(rep(1,length(regulParam)))
		    loss <- loss + vars[,k] + t(rep(1,p)) %*% bias[,,k]^2 
  		}

			#Estimation des directions EDR avec le paramètre de régularisation optimal
			regulOpt <- which(loss==min(loss))
			temp <- .edrNorm(M, Sigma + regulParam[regulOpt] * diag(p), K)
			matEDR <- temp$matEDR
			eigvalEDR <- temp$eigvalEDR
			additionalContent <- list(s=regulParam[regulOpt], estMSE=loss)

		}	
		if (method == "SR-SIR") {
			temp <- sliceMat(Y, X, H, details=TRUE)
			Xh <- t(temp$Xh)
			Ph <- temp$Ph

			testedEDR <- list()
			iterVec <- NULL
			rssVec <- NULL
			nbParamVec <- NULL
			gcvVec <- NULL							
			for (i in 1:length(regulParam)) {
				iter <- 0
				testedEDR[[i]] <- initEDR
				GNew <- Inf
				GDiff <- Inf
				
				#Alternating least squares algorithm				
				while(all(iter < maxIter, GDiff > tol)) {
					C1 <- solve(t(testedEDR[[i]]) %*% Sigma %*% Sigma %*% testedEDR[[i]]) %*% 
						t(testedEDR[[i]]) %*% Sigma
					C <- C1 %*% Xh
					C2 <- solve(C %*% Ph %*% t(C) %x% (Sigma %*% Sigma) + 
						regulParam[i] * diag(p * K))
					C3 <- C2 %*% ((C %*% Ph) %x% Sigma)
					testedEDR[[i]] <- matrix(C3 %*% as.vector(Xh), nrow=p, byrow=FALSE)
					
					S <- ((sqrt(Ph) %*% t(C)) %x% Sigma) %*% C2 %*% ((C %*% sqrt(Ph)) %x% Sigma)
					C4 <- (diag(p * H) - S) %*% (sqrt(Ph) %x% diag(p)) %*% as.vector(Xh)
					GOld <- GNew					
					GNew <- as.vector(t(C4) %*% C4)
					GDiff <- GOld - GNew					
					iter <- iter + 1
				}
				testedEDR[[i]] <- testedEDR[[i]] %*% 
					sqrt(diag(1/diag(t(testedEDR[[i]]) %*% testedEDR[[i]]), nrow=K, ncol=K))

				#Calcul des critères de sélection de paramètres
				iterVec <- c(iterVec, iter)				
				rssVec <- c(rssVec, GNew)
				pLambda <- sum(diag(S))
				nbParamVec <- c(nbParamVec, pLambda)
				gcvVec <- c(gcvVec, GNew / (p * H * (1 - pLambda / p / H) ^ 2))
			}
			chosenPar <- which(gcvVec == min(gcvVec))[1]
			matEDR <- testedEDR[[chosenPar]]
			eigvalEDR <- NULL
			additionalContent <- list(s=regulParam[chosenPar], testedEDR=testedEDR,
				iter=iterVec, gcv=gcvVec)
		}	

		res <- c(list(matEDR=matEDR, eigvalEDR=eigvalEDR, K=K, H=H, n=n, 
								method=paste("SIR-I,",method) , X=X, Y=Y), additionalContent)
	}	

	class(res) <-"edr"  
	res	

}


