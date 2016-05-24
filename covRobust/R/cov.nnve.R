"cov.nnve" <-
function(datamat, k = 12, pnoise = 0.050000000000000003, emconv = 0.001, bound = 1.5, extension = TRUE, devsm = 0.01)
{
	#
        #
        # Function to perform Nearest Neighbor Variance Estimation
        # Wang and Raftery(2002), "Nearest neighbor variance estimation (NNVE):
        #     Robust covariance estimation via nearest neighbor cleaning
        #                 (with discussion)", 
        #     Journal of the American Statistical Association 97:994-1019
        # Available as Technical Report 368 (2000) from
        #    http://www.stat.washington.edu/www/research/reports
	#
	# The data need to be a matrix of points in the process.
	# Each column contains one variable.
	#
	#
	#  Returns a list with the following components:
	#      cov            - robust covariance estimate
	#      mu             - mean
	#      postprob       - posterior probability
	#      classification - classification (0=noise otherwise 1)
	#                      (obtained by rounding postprob)
	#      innc - list of initial nearest-neighbor results 
        #            (components are the same as above)
	#
	#       
	datamat <- as.matrix(datamat)
	d <- dim(datamat)[2]
	n <- dim(datamat)[1]
        #
	#  The Re-scale NNC function
	#
	nclean.sub <- function(datamat, k, distances = NULL, convergence = 0.001,
		S.mean = NULL, S.sd = NULL)
	{
		#
		# Function to perform the Nearest Neighbour cleaning of
		#--------------------------------------------------------------
		# find the density of D_k
		#
		dDk <- function(x, lambda, k, d, alpha.d)
		{
			(exp( - lambda * alpha.d * x^d + log(2) + k * log(
				lambda * alpha.d) + log(x) * (d * k - 1) - log(
				gamma(k))))
		}
		# end of dDk
		#--------------------------------------------------------------
		# nearest-neighbor in S-PLUS
		#
		splusNN <- function(datamat, k)
		{
			n <- nrow(datamat)
			distances <- dist(datamat)
			#
			#  This next part sorts through the Splus distance object and forms kNNd, 
			#  kth nearest neighbour distance, for each point.
			#
			kNNd <- rep(0, n)
			N <- (n - 1):0
			I <- c(0, cumsum(N[-1]))
			J <- c(0, I + n - 1)
			a <- z <- NULL
			for(j in 1:n) {
				if(j > 1)
					a <- i + I[1:i]
				if(j < n)
					z <- J[j] + 1:N[j]
				kNNd[j] <- sort(distances[c(a, z)])[k]
				i <- j
			}
			kNNd
		}
		# end of splusNN
		#--------------------------------------------------------------
		d <- dim(datamat)[2]
		n <- dim(datamat)[1]
		#
		kthNND <- splusNN(t((t(datamat) - S.mean)/S.sd), k = k)
		alpha.d <- (2 * pi^(d/2))/(d * gamma(d/2))
		#
		# Now use kthNND in E-M algorithm, first get starting guesses.
		#
		delta <- rep(0, n)
		delta[kthNND > (min(kthNND) + diff(range(kthNND))/3)] <- 1
		p <- 0.5
		lambda1 <- k/(alpha.d * mean((kthNND[delta == 0])^d))
		lambda2 <- k/(alpha.d * mean((kthNND[delta == 1])^d))
		loglik.old <- 0
		loglik.new <- 1
		#
		# Iterator starts here, 
		#
		while(abs(loglik.new - loglik.old)/(1+abs(loglik.new)) > convergence) 
			{
			# E - step
			delta <- (p * dDk(kthNND, lambda1, k = k, d = d, 
				alpha.d = alpha.d))/(p * dDk(kthNND, lambda1,
				k = k, d = d, alpha.d = alpha.d) + (1 - p) *
				dDk(kthNND, lambda2, k = k, d = d, alpha.d = 
				alpha.d))
			# M - step
			p <- sum(delta)/n
			lambda1 <- (k * sum(delta))/(alpha.d * sum((kthNND^
				d) * delta))
			lambda2 <- (k * sum((1 - delta)))/(alpha.d * sum((
				kthNND^d) * (1 - delta)))
			loglik.old <- loglik.new
			loglik.new <- sum( - p * lambda1 * alpha.d * ((kthNND^
				d) * delta) - (1 - p) * lambda2 * alpha.d *
				((kthNND^d) * (1 - delta)) + delta * k * log(
				lambda1 * alpha.d) + (1 - delta) * k * log(
				lambda2 * alpha.d))
		}
		#
		# z will be the classifications. 1= in cluster. 0= in noise.
		#
		probs <- dDk(kthNND, lambda1, k = k, d = d, alpha.d = alpha.d)/
			(dDk(kthNND, lambda1, k = k, d = d, alpha.d = alpha.d) +
			dDk(kthNND, lambda2, k = k, d = d, alpha.d = alpha.d))
		mprob <- 1. - probs
		mu1 <- apply((probs * datamat), 2, sum)/sum(probs)
		mu2 <- apply((mprob * datamat), 2, sum)/sum(mprob)
		tpsig1 <- t(datamat) - mu1
		tpsig2 <- t(datamat) - mu2
		Sig1 <- tpsig1 %*% (probs * t(tpsig1))/sum(probs)
		Sig2 <- tpsig2 %*% (mprob * t(tpsig2))/sum(mprob)
		sd1 <- sqrt(diag(Sig1))
		sd2 <- sqrt(diag(Sig2))
		ans <- rbind(mu1, sd1, mu2, sd2)
		zz <- list(z = round(probs), kthNND = kthNND, probs = probs,
			p = p, mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2,
			lambda1 = lambda1, lambda2 = lambda2, Sig1 = Sig1,
			Sig2 = Sig2, ans = ans)
		return(zz)
	}
	#----- end of nclean.sub-----------------------------------------------
	#
	#    Two procedures to link between a symmetric matrix and its vec(.)
	#
	#
	Mtovec <- function(M)
	{
		n <- dim(M)[1]
		d <- dim(M)[2]
		if(abs(n - d) > 0.01) {
			cat("The input has to be a square matrix")
		}
		else {
			vec <- rep(0, 0)
			for(i in (1:n)) {
				for(j in (i:d)) {
					vec <- c(vec, M[i, j])
				}
			}
			vec
		}
	}
	vectoM <- function(vec, d)
	{
		n <- length(vec)
		M <- matrix(rep(0, d * d), d, d)
		L <- 1
		for(i in 1:d) {
			for(j in i:d) {
				M[i, j] <- vec[L]
				L <- L + 1
				M[j, i] <- M[i, j]
			}
		}
		M
	}
	#
	#----------------------------------------------------------------------
	#
	pd <- dim(datamat)[2]
	S.mean <- apply(datamat, 2, median)
	S.sd <- apply(datamat, 2, mad)
	#
	#  NNC based on original data
	#
	orgNNC <- nclean.sub(datamat, k, convergence = 0.001, S.mean = S.mean,
		S.sd = S.sd)
	nnoise <- min(c(sum(1 - orgNNC$z), round(pnoise * n)))
	knnd <- orgNNC$kthNND
	ord <- (n + 1) - rank(knnd)
	muT <- orgNNC$mu1
	SigT <- orgNNC$Sig1
	SigT <- (SigT + t(SigT))/2.
	SigTN <- diag(orgNNC$sd1^2)
	if(nnoise > 6) {
		ncho <- nnoise
		ncho1 <- floor(ncho/2)
		ncho2 <- ncho - ncho1
		cho <- (1:n)[ord <= ncho1]
		xcho <- datamat[cho,  ]
		ev <- eigen(SigT)
		evv <- ev$values
		minv <- max((1:d)[evv > 9.9999999999999998e-13])
		if(minv > 2) {
			vv1 <- ev$vectors[, (minv - 1)]
			vv2 <- ev$vectors[, minv]
		}
		else {
			vv1 <- ev$vectors[, 1]
			vv2 <- ev$vectors[, 2]
		}
		ot <- acos(sum(vv1 * vv2)/(sum(vv1^2) * sum(vv2^2))^0.5)
		for(kk1 in 1:(ncho2)) {
			pseg <- 1/(ncho2 + 1) * kk1 * ot
			xcho <- rbind(xcho, (sin(pseg) * vv1 + cos(pseg) * vv2 +
				muT))
		}
	}
	else {
		nnoise <- 3
		cho <- (1:n)[ord <= nnoise]
		xcho <- datamat[cho,  ]
	}
	n2 <- (dim(xcho))[1]
	schox <- mahalanobis(xcho, muT, SigTN)
	Nc <- matrix(rep(muT, n2), nrow = n2, byrow = TRUE)
	Ndir <- (xcho - Nc)/(schox^0.5)
	#
	# initial set up
	#
	ch1 <- c(qchisq(0.98999999999999999, pd), qchisq(1 - 10^(-4), pd))
	Xa <- seq(ch1[1], ch1[2], length = 6)
	gap <- Xa[2] - Xa[1]
	initv <- diag(orgNNC$Sig1)
	xa <- Xa[1]
	SaveM <- c(xa, orgNNC$mu1, Mtovec(orgNNC$Sig1))
	OldP <- orgNNC$probs
	SaveP <- OldP
	Np <- Nc - Ndir * (xa^0.5)
	updNNC <- nclean.sub(rbind(datamat, Np), k, convergence = 0.001, S.mean = S.mean,
		S.sd = S.sd)
	SaveM <- rbind(SaveM, c(xa, updNNC$mu1, Mtovec(updNNC$Sig1)))
	SaveP <- rbind(SaveP, (updNNC$probs)[1:n])
	#sda<-Mtovec(orgNNC$Sig1)  
	#sda save the results corresponding to xa = qchisq(.99,pd)
	stopv <- diag(updNNC$Sig1)
	time1 <- 2
	#
	#
	while((time1 <= 6) && (all(stopv < (1 + bound) * initv))) {
		xa <- Xa[time1]
		Np <- Nc - Ndir * (xa^0.5)
		updNNC <- nclean.sub(rbind(datamat, Np), k, convergence = 0.001, S.mean = 
			S.mean, S.sd = S.sd)
		SaveM <- rbind(SaveM, c(xa, updNNC$mu1, Mtovec(updNNC$Sig1)))
		SaveP <- rbind(SaveP[2,  ], (updNNC$probs)[1:n])
		time1 <- time1 + 1
		stopv <- diag(updNNC$Sig1)
		NULL
	}
	#
	#
	#   Procedure stop if the added noise cause a "surge" within the range
	#   sdb save the results within the given "range"
	#
	#
	if(all(stopv < (1 + bound) * initv)) {
		dSaveM <- dim(SaveM)[1]
		ans <- SaveM[dSaveM,  ]
		sdb <- SaveM[dSaveM,  ]
		NewP <- SaveP[2,  ]
		#
		#  adding extension
		#
		if(extension) {
			time2 <- 1
			Fstop <- FALSE
			tpv <- stopv
			while((time2 < 2) && (all(stopv < (1 + bound) * initv))
				) {
				xa <- xa + gap
				startv <- stopv
				Np <- Nc - Ndir * (xa^0.5)
				updNNC <- nclean.sub(rbind(datamat, Np), k, convergence = 
					0.001, S.mean = S.mean, S.sd = S.sd)
				SaveM <- rbind(SaveM, c(xa, updNNC$mu1, Mtovec(
					updNNC$Sig1)))
				SaveP <- rbind(SaveP[2,  ], (updNNC$probs)[
					1:n])
				stopv <- apply(rbind((startv * 2 - tpv), diag(
					updNNC$Sig1)), 2, mean)
				tpv <- diag(updNNC$Sig1)
				Fstop <- all((abs(stopv - startv) <= ((1+abs(startv)) *
					devsm)))
				if(Fstop)
					time2 <- time2 + 1
				else time2 <- 1
				NULL
			}
			#end-while
			# checking the stop criterior at the end of extension
			if(all(stopv < (1 + bound) * initv)) {
				dSaveM <- dim(SaveM)[1]
				ans <- SaveM[dSaveM,  ]
				NewP <- SaveP[2,  ]
			}
			else {
				dSaveM <- dim(SaveM)[1]
				ans <- SaveM[dSaveM - 1,  ]
				NewP <- SaveP[1,  ]
			}
		}
	}
	else {
		dSaveM <- dim(SaveM)[1]
		ans <- SaveM[dSaveM - 1,  ]
		sdb <- ans[-1]
		NewP <- SaveP[1,  ]
	}
	nncvar <- vectoM(ans[ - (1:(1 + pd))], pd)
	mu <- ans[2:(1 + pd)]
	#       
	list(cov = nncvar, mu = mu, postprob = NewP, classification = round(NewP), innc = list(cov = orgNNC$Sig1, mu = orgNNC$mu1, postprob = OldP, classification = round(OldP)))
}
