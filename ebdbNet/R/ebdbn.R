`ebdbn` <-
function(y, K, input = "feedback", conv.1 = .15, conv.2 = .05, conv.3 = .01, verbose = TRUE, max.iter = 100, max.subiter = 200)
{

	## M is input dimension (= P if feedback)
	## alpha and gamma are of dimension K, beta and delta are of dimension M, v is of dimension P

	if(is.list(y) != TRUE) {
      	stop("Error: ", paste(sQuote("y"), sep = ""), " should be a list.")
    	}

	R <- length(y)
	P <- dim(y[[1]])[1]
	T <- dim(y[[1]])[2]

	if(is.list(input) == FALSE) {
		u <- fdbkFunc(y)
		M <- P
	}
	if(is.list(input) == TRUE) {
		u <- input
		M <- dim(u[[1]])[1]
		if(length(u) != R) {
			stop("Error: ", paste(sQuote("input"), sep = ""), " should be a list of length R of MxT matrices for input networks.")
		}
		if(dim(u[[1]])[2] != T) {
			stop("Error: ", paste(sQuote("input"), sep = ""), " should be a list of length R of MxT matrices for input networks.")
		}
	}

	## Initial values of x.0, alpha, beta, gamma, delta, v, mu, sigma
	if(K > 0) {
		x.0 <- vector("list", R)
		for(r in 1:R) {x.0[[r]] <- matrix(rnorm(K*T), nrow = K, ncol = T)}
		alpha.0 <- runif(K)
		beta.0 <- runif(M)
		gamma.0 <- runif(K)
		mu.0 <- rep(0, K)
		sigma.0 <- rep(1, K)
	}
	if(K == 0) {
		x.0 <- alpha.0 <- beta.0 <- gamma.0 <- mu.0 <- sigma.0 <- 0
	}
	delta.0 <- runif(M)
	v.0 <- runif(P)


	APost <- rep(0, K*K)
	BPost <- rep(0, M*K)
	CPost <- rep(0, P*K)
	CvarPost <- rep(0, P*K*K)
	DPost <- rep(0, P*M)
	DvarPost <- rep(0, M*M*P)

	if(K>0) {x0 <- as.vector(t(do.call(rbind, x.0)))}
	if(K==0) {x0 <- 0}
	yorig <- as.vector(t(do.call(rbind, y)))
	uorig <- as.vector(t(do.call(rbind, u)))

	verboseInd <- ifelse(verbose == TRUE, 1, 0)

	## Run Full Algorithm and Find Posterior Mean and Variance
	test = .C("RunWrapGen", R = as.integer(R), P = as.integer(P), 
		T = as.integer(T), K = as.integer(K), M = as.integer(M), xx = as.double(x0), 
		yy = as.double(yorig), uu = as.double(uorig), alpha = as.double(alpha.0), 
		beta = as.double(beta.0), gamma = as.double(gamma.0), 
		delta = as.double(delta.0), v = as.double(v.0),
		mu = as.double(mu.0), sigma = as.double(sigma.0),
		conv1 = as.double(conv.1), conv2 = as.double(conv.2),
		conv3 = as.double(conv.3), APost = as.double(APost),
		BPost = as.double(BPost), CPost = as.double(CPost),
		DPost = as.double(DPost), CvarPost = as.double(CvarPost),
		DvarPost = as.double(DvarPost), alliterations = as.integer(0), 
		maxiterations = as.integer(max.iter), subiterations = as.integer(max.subiter), 
		verboseInd = as.integer(verboseInd), PACKAGE = "ebdbNet")

	## Estimated Posterior Mean and Variance

	APost <- matrix(test$APost, nrow = K, ncol = K, byrow = TRUE)
	BPost <- matrix(test$BPost, nrow = K, ncol = M, byrow = TRUE)
	CPost <- matrix(test$CPost, nrow = P, ncol = K, byrow = TRUE)
	DPost <- matrix(test$DPost, nrow = P, ncol = M, byrow = TRUE)
	CvarPost <- vector("list", P)
	DvarPost <- vector("list", P)
	for(i in 1:P) {
		DvarPost[[i]] <- matrix(test$DvarPost[(M*M*(i-1)+1):(M*M*i)], nrow = M, ncol = M, byrow = TRUE)
		CvarPost[[i]] <- matrix(test$CvarPost[(K*K*(i-1)+1):(K*K*i)], nrow = K, ncol = K, byrow = TRUE)
	}

	## Estimated hidden states
	xPost <- vector("list", R)
	for(r in 1:R) {
		xPost[[r]] <- matrix(test$xx[(K*T*(r-1)+1):(K*T*r)], nrow = K, ncol = T, byrow = TRUE)
	}

	## Estimated hyperparameters
	alphaEst <- test$alpha
	betaEst <- test$beta
	gammaEst <- test$gamma
	deltaEst <- test$delta
	vEst <- test$v
	muEst <- test$mu
	sigmaEst <- test$sigma

	## Total number of iterations needed
	alliterations <- test$alliterations

	## Calculate z-statistics for D
	Dvar <- matrix(0, nrow = P, ncol = M)
	for(i in 1:P) Dvar[i,] <- diag(DvarPost[[i]]);
	Dsd <- sqrt(Dvar)
	z <- DPost / Dsd

	rownames(DPost) <- rownames(z) <- rownames(y[[1]])
	if(is.list(input) == TRUE) {
		colnames(DPost) <- colnames(z) <- rownames(u[[1]])
	}
	if(is.list(input) == FALSE) {
		colnames(DPost) <- colnames(z) <- rownames(y[[1]])
	}
	type <- ifelse(is.list(input) == TRUE, "input", "feedback")


	results <- list(APost = APost, BPost = BPost, CPost = CPost, DPost = DPost, CvarPost = CvarPost, 
		DvarPost = DvarPost, xPost = xPost, alphaEst = alphaEst, betaEst = betaEst, 
		gammaEst = gammaEst, deltaEst = deltaEst, vEst = vEst, muEst = muEst, sigmaEst = sigmaEst,
		alliterations = alliterations, z = z, type = type)

	class(results) <- "ebdbNet"
	return(results)
}

