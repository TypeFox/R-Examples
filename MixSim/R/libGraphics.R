
RotatedMixture <- function(Pi, Mu, S){

	K <- dim(Mu)[1]
	p <- dim(Mu)[2]

	SigmaMix <- matrix(rep(0, p * p), ncol = p)
	WeiMuVec <- rep(0, p)

	RotMu <- Mu
	RotS <- S
	

	for (i in 1:K){

	    SigmaMix <- SigmaMix + Pi[i] * (S[,,i] + Mu[i,] %*% t(Mu[i,]))

	    WeiMuVec <- WeiMuVec + Pi[i] * Mu[i,]

	}

	WeiMuMat <- WeiMuVec %*% t(WeiMuVec)
 
	SigmaMix <- SigmaMix - WeiMuMat
	
	
	Gamm <- eigen(SigmaMix)$vectors

	for (i in 1:K){

	    RotMu[i,] <- t(Gamm) %*% Mu[i,]
	    RotS[,,i] <- t(Gamm) %*% S[,,i] %*% Gamm

	}


	return(list(RotMu = RotMu, RotS = RotS))


}


RotatedStandMixture <- function(Pi, Mu, S){

	K <- dim(Mu)[1]
	p <- dim(Mu)[2]

	SigmaMix <- matrix(rep(0, p * p), ncol = p)
	WeiMuVec <- rep(0, p)

	RotMu <- Mu
	RotS <- S
	

	for (i in 1:K){

	    SigmaMix <- SigmaMix + Pi[i] * (S[,,i] + Mu[i,] %*% t(Mu[i,]))

	    WeiMuVec <- WeiMuVec + Pi[i] * Mu[i,]

	}

	WeiMuMat <- WeiMuVec %*% t(WeiMuVec)
 
	SigmaMix <- SigmaMix - WeiMuMat
	
	Eig <- eigen(SigmaMix)
	
	Gamm <- Eig$vectors
	LambdaNH <- diag(1/sqrt(Eig$values))

	for (i in 1:K){

	    RotMu[i,] <- LambdaNH %*% t(Gamm) %*% Mu[i,]
	    RotS[,,i] <- LambdaNH %*% t(Gamm) %*% S[,,i] %*% Gamm %*% LambdaNH

	}


	return(list(RotMu = RotMu, RotS = RotS, SigmaMix = SigmaMix, Lambda = diag(eigen(SigmaMix)$values)))


}



MultiPolygons <- function(xx, yy, colors, Nx){

	xx1 <- xx[1]
	xx2 <- xx[2]
	yy11 <- yy[1]
	yy12 <- yy[2]
	yy21 <- yy[3]
	yy22 <- yy[4]

	stepx <- (xx2 - xx1) / Nx
	stepy1 <- (yy12 - yy11) / Nx
	stepy2 <- (yy22 - yy21) / Nx

	x0 <- xx1
	y10 <- yy11
	y20 <- yy21

	for (i in 1:Nx){

		x1 <- x0 + stepx
		y11 <- y10 + stepy1
		y21 <- y20 + stepy2

		polygon(c(x0, x1, x1, x0), c(y20, y21, y11, y10), col = colors[i], border = NA)

		x0 <- x1
		y10 <- y11
		y20 <- y21

	}

}


ConstractPDPlot <- function(x, Pi, sd, file, Ny, Nx, MaxInt, marg){

	x <- t(x)

	re <- c(1, 0, 0, 1, 0, 1, 0.3, 1, 1, 0.7, 0.3, 0.7, 0.3, 0, 0, 0.7, 0.3, 0.7)
	gr <- c(0, 1, 0, 0, 1, 1, 0.7, 0.7, 0.3, 0.3, 1, 1, 0.7, 0.7, 0.3, 0.3, 0, 0)
	bl <- c(0, 0, 1, 1, 1, 0, 1, 0.3, 0.7, 1, 0.7, 0.3, 0, 0.3, 0.7, 0, 0.7, 0.3)

	p <- dim(x)[1]
	K <- dim(x)[2]

	hei <- array(rep(NA, Ny * K * p), c(Ny, p, K))
	tr.lev <- array(rep(NA, Ny * K * (p-1) * Nx), c(Ny, (p-1) * Nx, K))

	Z <- seq(0, 1, length.out = Ny + 2)[-1]
	Z <- qnorm(1 - Z / 2)[(Ny+1):1]

	lZ <- Z[-(Ny + 1)]
	uZ <- Z[-1]

	mZ <- (lZ + uZ) / 2

	for (j in 1:K){
		for (h in 1:p){
			for (i in 1:Ny){
				hei[i,h,j] <- dnorm(mZ[i] * sd[h,j], sd = sd[h,j])
			}
		}
	}

	hei <- MaxInt / max(hei) * hei

	for (j in 1:K){
		for (h in 1:(p - 1)){
			for (i in 1:Ny){
				tr.lev[i,((h-1)*Nx+1):(h*Nx),j] <- seq(hei[i,h,j], hei[i,h+1,j], length.out = Nx+2)[c(-1,-(Nx+2))]
			}
		}
	}


	color <- array(rep(NA, K * Ny * (p-1)*Nx), c(Ny, (p-1)*Nx, K))

	PiMax <- max(Pi)

	coll <- rep(NA,K)

	for (i in 1:Ny){
		for (g in 1:(Nx*(p-1))){
			for (k in 1:K){
			    color[i,g,k] <- rgb(re[k], gr[k], bl[k], tr.lev[i,g,k] * Pi[k] / PiMax)
			}	
		}
	}


	for (k in 1:K){
		coll[k] <- rgb(re[k], gr[k], bl[k], 1)
	}

	if (!is.null(file)){
		pdf(file = file, version = "1.4")
	}
	par(mar = c(marg[1], marg[2], marg[3], marg[4]))

	plot( c(1, p), c( min(x - max(uZ) * sd), max(x + max(uZ) * sd) ), type = "n", xlab
	= "", ylab = "", axes = FALSE)

	axis(1, tick = TRUE, xaxp = c(1, p, p - 1))
	box()


	for (j in 1:K){

	    lines(1:p, x[,j], col = coll[j], lty = 1, lwd = .1)   
	    lines(1:p, x[,j], col = 1, lty = 2, lwd = .1) 

	    for (i in 1:Ny){

		for (h in 1:(p-1)){

			yy11 <- x[h,j] - uZ[i] * sd[h,j]
			yy12 <- x[h+1,j] - uZ[i] * sd[h+1,j]
			yy21 <- x[h,j] - lZ[i] * sd[h,j]
			yy22 <- x[h+1,j] - lZ[i] * sd[h+1,j]

			MultiPolygons(c(h, h+1), c(yy11, yy12, yy21, yy22), color[i,((h-1)*Nx+1):(h*Nx),j], Nx)


			yy11 <- x[h,j] + uZ[i] * sd[h,j]
			yy12 <- x[h+1,j] + uZ[i] * sd[h+1,j]
			yy21 <- x[h,j] + lZ[i] * sd[h,j]
			yy22 <- x[h+1,j] + lZ[i] * sd[h+1,j]

			MultiPolygons(c(h,h+1), c(yy11, yy12, yy21, yy22), color[i,((h-1)*Nx+1):(h*Nx),j], Nx)

		}

	    }

	}

	if (!is.null(file)){
		dev.off()
	}

}

pdplot <- function(Pi, Mu, S, file = NULL, Nx = 5, Ny = 5, MaxInt = 1.0, marg = c(2,1,1,1)){

	if (sum((Pi <= 0) | (Pi >= 1)) != 0) stop("Wrong vector of mixing proportions Pi...\n")
	if (Nx < 1) stop("Wrong value of Nx...\n")
	if (Ny < 1) stop("Wrong value of Ny...\n")
	if ((MaxInt < 0) | (MaxInt > 1)) stop("Wrong value of MaxInt...\n")

	K <- dim(Mu)[1]
	p <- dim(Mu)[2]

	RM <- RotatedMixture(Pi, Mu, S)

	sd <- array( rep(NA, p * K), c(p, K) )
	for (k in 1:K){
		sd[,k] <- sqrt(diag(RM$RotS[,,k]))
	}

	ConstractPDPlot(RM$RotMu, Pi, sd, file, Ny, Nx, MaxInt, marg)

}

