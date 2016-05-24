
overlap <- function(Pi, Mu, S, eps = 1e-06, lim = 1e06){
	
	if (sum((Pi <= 0) | (Pi >= 1)) != 0) stop("Wrong vector of mixing proportions Pi...\n")
	if (eps <= 0) stop("Wrong value of eps...\n")
	if (lim < 1) stop("Wrong value of lim...\n")

	K <- dim(Mu)[1]
	p <- dim(Mu)[2]

        Mu1 <- as.vector(t(Mu))
        S1 <- as.vector(S)
        OmegaMap1 <- rep(0, K*K)

        rcMax <- c(0, 0)
        

	Q <- .C("runExactOverlap", p1 = as.integer(p), K1 = as.integer(K), Pi = as.double(Pi), Mu1 = as.double(Mu1), S1 = as.double(S1), pars = as.double(c(eps, eps)), lim1 = as.integer(lim), OmegaMap1 = as.double(OmegaMap1), BarOmega = as.double(1), MaxOmega = as.double(1), EigOmega = as.double(1), rcMax = as.integer(rcMax), PACKAGE = "MixSim")

        return(list(OmegaMap = matrix(Q$OmegaMap1, byrow = TRUE, ncol = K), BarOmega = Q$BarOmega, MaxOmega = Q$MaxOmega, rcMax = Q$rcMax + 1))

}



overlapGOM <- function(Pi, Mu, S, eps = 1e-06, lim = 1e06){
	
	if (sum((Pi <= 0) | (Pi >= 1)) != 0) stop("Wrong vector of mixing proportions Pi...\n")
	if (eps <= 0) stop("Wrong value of eps...\n")
	if (lim < 1) stop("Wrong value of lim...\n")

	K <- dim(Mu)[1]
	p <- dim(Mu)[2]

        Mu1 <- as.vector(t(Mu))
        S1 <- as.vector(S)
        OmegaMap1 <- rep(0, K*K)

        rcMax <- c(0, 0)
        

	Q <- .C("runExactOverlap", p1 = as.integer(p), K1 = as.integer(K), Pi = as.double(Pi), Mu1 = as.double(Mu1), S1 = as.double(S1), pars = as.double(c(eps, eps)), lim1 = as.integer(lim), OmegaMap1 = as.double(OmegaMap1), BarOmega = as.double(1), MaxOmega = as.double(1), EigOmega = as.double(1), rcMax = as.integer(rcMax), PACKAGE = "MixSim")

        return(Q$EigOmega)

}



MixSim <- function(BarOmega = NULL, MaxOmega = NULL, K, p, sph = FALSE, hom = FALSE, ecc = 0.90, PiLow = 1.0, int = c(0.0, 1.0), resN = 100, eps = 1e-06, lim = 1e06){

	if (p < 1) stop("Wrong number of dimensions p...\n")
	if (K <= 1) stop("Wrong number of mixture components K...\n")
	if ((sph != FALSE) & (sph != TRUE)) stop("Wrong value of sph p...\n")
	if ((hom != FALSE) & (hom != TRUE)) stop("Wrong value of hom p...\n")
	if ((ecc <= 0) | (ecc > 1)) stop("Wrong value of ecc...\n")
	if ((PiLow <= 0) | (PiLow > 1)) stop("Wrong value of PiLow...\n")
	if (int[1] >= int[2]) stop("Wrong interval int...\n")
	if (resN < 1) stop("Wrong value of resN...\n")
	if (eps <= 0) stop("Wrong value of eps...\n")
	if (lim < 1) stop("Wrong value of lim...\n")


        Pi <- rep(0, K)
        Mu1 <- rep(0, K*p)
        S1 <- rep(0, K*p*p)
        OmegaMap1 <- rep(0, K*K)
        rcMax <- c(0, 0)
	Lbound <- int[1]
	Ubound <- int[2]

 
        if ((is.null(MaxOmega)) & (!is.null(BarOmega))){
               method <- 0
               Omega <- BarOmega
        }
        if ((is.null(BarOmega)) & (!is.null(MaxOmega))){
               method <- 1
               Omega <- MaxOmega
             }
        if ((!is.null(BarOmega)) & (!is.null(MaxOmega)))  method <- 3   
	if ((is.null(BarOmega)) & (is.null(MaxOmega)))  method <- -1   

        if ((method == 0) | (method == 1)){

              Q <- .C("runOmegaClust", Omega1 = as.double(Omega), method1 = as.integer(method), p1 = as.integer(p), K1 = as.integer(K), PiLow1 = as.double(PiLow), Lbound1 = as.double(Lbound), Ubound1 = as.double(Ubound), emax1 = as.double(ecc), pars = as.double(c(eps, eps)), lim1 = as.integer(lim), resN1 = as.integer(resN), sph1 = as.integer(sph), hom1 = as.integer(hom),Pi = as.double(Pi), Mu1 = as.double(Mu1), S1 = as.double(S1), OmegaMap1 = as.double(OmegaMap1), BarOmega = as.double(1), MaxOmega = as.double(1), EigOmega = as.double(1), rcMax = as.integer(rcMax), fail = as.integer(1), PACKAGE = "MixSim")


        }

        if (method == 3){
              

              Q <- .C("runOmegaBarOmegaMax", p1 = as.integer(p), K1 = as.integer(K), PiLow1 = as.double(PiLow), Lbound1 = as.double(Lbound), Ubound1 = as.double(Ubound), emax1 = as.double(ecc), pars = as.double(c(eps, eps)), lim1 = as.integer(lim), resN1 = as.integer(resN), sph1 = as.integer(sph), Pi = as.double(Pi), Mu1 = as.double(Mu1), S1 = as.double(S1), OmegaMap1 = as.double(OmegaMap1), BarOmega = as.double(BarOmega), MaxOmega = as.double(MaxOmega), rcMax = as.integer(rcMax), fail = as.integer(1), PACKAGE = "MixSim")

        }

        if (method != -1){

		if (Q$fail == 0){
			ret <- list(Pi = Q$Pi, Mu = matrix(Q$Mu1, byrow = TRUE, ncol = p), S = array(Q$S1, c(p, p, K)), OmegaMap = matrix(Q$OmegaMap1, byrow = TRUE, ncol = K), BarOmega = Q$BarOmega, MaxOmega = Q$MaxOmega, rcMax = Q$rcMax + 1, fail = Q$fail)
			class(ret) <- "MixSim"
			return(ret)
		}

	} else {
		cat("Error: at least one overlap characteristic should be specified...\n")		
	}
		
}


MixGOM <- function(goMega = NULL, K, p, sph = FALSE, hom = FALSE, ecc = 0.90, PiLow = 1.0, int = c(0.0, 1.0), resN = 100, eps = 1e-06, lim = 1e06){

	if (p < 1) stop("Wrong number of dimensions p...\n")
	if (K <= 1) stop("Wrong number of mixture components K...\n")
	if ((sph != FALSE) & (sph != TRUE)) stop("Wrong value of sph p...\n")
	if ((hom != FALSE) & (hom != TRUE)) stop("Wrong value of hom p...\n")
	if ((ecc <= 0) | (ecc > 1)) stop("Wrong value of ecc...\n")
	if ((PiLow <= 0) | (PiLow > 1)) stop("Wrong value of PiLow...\n")
	if (int[1] >= int[2]) stop("Wrong interval int...\n")
	if (resN < 1) stop("Wrong value of resN...\n")
	if (eps <= 0) stop("Wrong value of eps...\n")
	if (lim < 1) stop("Wrong value of lim...\n")


        Pi <- rep(0, K)
        Mu1 <- rep(0, K*p)
        S1 <- rep(0, K*p*p)
        OmegaMap1 <- rep(0, K*K)
        rcMax <- c(0, 0)
	Lbound <- int[1]
	Ubound <- int[2]


        if (!is.null(goMega)){

		method <- 2

	        Q <- .C("runOmegaClust", Omega1 = as.double(goMega), method1 = as.integer(method), p1 = as.integer(p), K1 = as.integer(K), PiLow1 = as.double(PiLow), Lbound1 = as.double(Lbound), Ubound1 = as.double(Ubound), emax1 = as.double(ecc), pars = as.double(c(eps, eps)), lim1 = as.integer(lim), resN1 = as.integer(resN), sph1 = as.integer(sph), hom1 = as.integer(hom),Pi = as.double(Pi), Mu1 = as.double(Mu1), S1 = as.double(S1), OmegaMap1 = as.double(OmegaMap1), BarOmega = as.double(1), MaxOmega = as.double(1), EigOmega = as.double(goMega), rcMax = as.integer(rcMax), fail = as.integer(1), PACKAGE = "MixSim")

		if (Q$fail == 0){
			ret <- list(Pi = Q$Pi, Mu = matrix(Q$Mu1, byrow = TRUE, ncol = p), S = array(Q$S1, c(p, p, K)), goMega = Q$EigOmega, fail = Q$fail)
			class(ret) <- "MixGOM"
			return(ret)
		}

	} else {
		cat("Error: goMega value should be specified...\n")		
	}
		
}



### Add for revision.
print.MixSim <- function(x, ...){
        K <- length(x$Pi)
        p <- ncol(x$Mu)
	cat("K = ", K,
	    ", p = ", p,
            ", BarOmega = ", x$BarOmega,
            ", MaxOmega = ", x$MaxOmega,
            ", success = ", ifelse(x$fail == 0, "TRUE", "FALSE"),
            ".\n", sep = "")
        cat("\nPi: \n")
        print(x$Pi)
        cat("\nMu: \n")
        Mu <- x$Mu
        colnames(Mu) <- paste("p.", 1:p, sep = "")
        rownames(Mu) <- paste("K.", 1:K, sep = "")
        print(Mu)
        cat("\nS: ... too long and skipped. Use operator $ to access.\n")
	invisible()
} # End of print.MixSim().

summary.MixSim <- function(object, ...){
        K <- length(object$Pi)
        OmegaMap <- object$OmegaMap
        colnames(OmegaMap) <- paste("k.", 1:K, sep = "")
        rownames(OmegaMap) <- paste("k.", 1:K, sep = "")
	cat("OmegaMap: \n")
        print(OmegaMap)
	cat("\nrcMax:", object$rcMax, "\n")
	invisible()
} # End of summary.MixSim().

print.MixGOM <- function(x, ...){
        K <- length(x$Pi)
        p <- ncol(x$Mu)
	cat("K = ", K,
	    ", p = ", p,
            ", goMega = ", x$goMega,
            ", success = ", ifelse(x$fail == 0, "TRUE", "FALSE"),
            ".\n", sep = "")
        cat("\nPi: \n")
        print(x$Pi)
        cat("\nMu: \n")
        Mu <- x$Mu
        colnames(Mu) <- paste("p.", 1:p, sep = "")
        rownames(Mu) <- paste("K.", 1:K, sep = "")
        print(Mu)
        cat("\nS: ... too long and skipped. Use operator $ to access.\n")
	invisible()
} # End of print.MixGOM().



RandIndex <- function(id1, id2){

	if (length(id1) != length(id2)) stop("Lengths of partitioning vectors do not match...\n")

	n <- length(id1)
	

	A <- as.factor(id1)
	B <- as.factor(id2)
	K1 <- nlevels(A)
	K2 <- nlevels(B)

	for (i in 1:nlevels(A)){
		ind <- A == levels(A)[i]
		id1[ind] <- i - 1
	}
	for (i in 1:nlevels(B)){
		ind <- B == levels(B)[i]
		id2[ind] <- i - 1
	}


	Q <- .C("runAdjRand", n = as.integer(n), K1 = as.integer(K1), K2 = as.integer(K2), id1 = as.integer(id1), id2 = as.integer(id2), Rand = as.double(0), aRand = as.double(0), F = as.double(0), PACKAGE = "MixSim")

       	return(list(R = Q$Rand, AR = Q$aRand, F = Q$F, M = n * (n - 1) * (1 - Q$Rand)))

}





ClassProp <- function(id1, id2){

	if (length(id1) != length(id2)) stop("Lengths of partitioning vectors do not match...\n")

	n <- length(id1)
	

	A <- as.factor(id1)
	B <- as.factor(id2)
	K1 <- nlevels(A)
	K2 <- nlevels(B)

	for (i in 1:nlevels(A)){
		ind <- A == levels(A)[i]
		id1[ind] <- i - 1
	}
	for (i in 1:nlevels(B)){
		ind <- B == levels(B)[i]
		id2[ind] <- i - 1
	}

	if (min(K1, K2) == 1) return(max(table(id1, id2)))

	Q <- .C("runProAgree", n = as.integer(n), K1 = as.integer(K1), K2 = as.integer(K2), id1 = as.integer(id1), id2 = as.integer(id2), maxPro = as.double(0), PACKAGE = "MixSim")

       	return(Q$maxPro)

}





VarInf <- function(id1, id2){

	if (length(id1) != length(id2)) stop("Lengths of partitioning vectors do not match...\n")

	n <- length(id1)
	

	A <- as.factor(id1)
	B <- as.factor(id2)
	K1 <- nlevels(A)
	K2 <- nlevels(B)

	for (i in 1:nlevels(A)){
		ind <- A == levels(A)[i]
		id1[ind] <- i - 1
	}
	for (i in 1:nlevels(B)){
		ind <- B == levels(B)[i]
		id2[ind] <- i - 1
	}


	Q <- .C("runVarInf", n = as.integer(n), K1 = as.integer(K1), K2 = as.integer(K2), id1 = as.integer(id1), id2 = as.integer(id2), VI = as.double(0), PACKAGE = "MixSim")

       	return(Q$VI)

}





perms <- function(n){

	n <- floor(n)
	if (n <= 0) stop("Incorrect number of elements...\n")
	if (n == 1) return(1)

	permN <- factorial(n)
	perms1 <- rep(0, n * permN)

	Q <- .C("runPerms", n1 = as.integer(n), permN1 = as.integer(permN), perms = as.integer(perms1), PACKAGE = "MixSim")

       	return(matrix(Q$perms, ncol = n, byrow = T) + 1)

}









