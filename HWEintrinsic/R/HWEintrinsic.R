hwe.ibf.mc <- function(y, t, M = 10000, verbose = TRUE) {
	if (class(y) != "HWEdata") {
		stop("y argument not of class 'HWEdata'. Type '?HWEdata' for help.")
	}
	y.mat <- y@data.mat
	r <- y@size
	R <- r*(r + 1)/2
	y.vec <- y@data.vec
	n <- sum(y.vec, na.rm = TRUE)
	if ((t <= 0) | (t > n) | (is.null(t))) {
		stop("The training sample size t has to be an integer in between 1 and the sample size n.")
	}
	r.i <- rowSums(y.mat, na.rm = TRUE)
	c.i <- colSums(y.mat, na.rm = TRUE)
	theta.hat <- (y.vec + 1)/(n + R)
	const <- lfactorial(t) + lfactorial(t + R - 1) + lfactorial(2*n + r - 1) - lfactorial(2*t + r - 1) - lfactorial(n + t + R - 1) -
			(n - sum(diag(y.mat)))*log(2) - sum(lfactorial(r.i + c.i))
	bf.draws <- vector(length = M)
	for (k in 1:M) {
		x.vec <- rmultinom(1, t, theta.hat)
		x <- matrix(NA, r, r)
		for (i in 1:r) {
			x[i, 1:i] <- x.vec[((i - 1)*i/2 + 1):(i*(i + 1)/2)]
		}
		z.i <- rowSums(x, na.rm = TRUE)
		d.i <- colSums(x, na.rm = TRUE)
		bf.draws[k] <- const + (t - sum(diag(x)))*log(2) + sum(lfactorial(z.i+d.i)) - 2*sum(lfactorial(x.vec)) +
				sum(lfactorial(x.vec + y.vec)) - dmultinom(x = x.vec, prob = theta.hat, log = TRUE)
		if (verbose) {
			if (round(k/10000) == k/10000) print(paste(k, "/", M, sep = ""), quote = FALSE)
		}
	}
	bf.draws <- exp(bf.draws)
	bf_mc <- mean(bf.draws)
	npp <- 1/(1 + bf_mc)
	out <- new("HWEintr", bf = bf_mc, npp = npp, draws = bf.draws, data = y.mat)
	return(out)
}

hwe.ibf <- function(y, t) {
	if (class(y) != "HWEdata") {
		stop("y argument not of class 'HWEdata'. Type '?HWEdata' for help.")
	}
	r <- y@size
	if (r != 2) {
		stop("The exact calculation for the Bayes Factor based on intrinsic priors is available only for r (number of alleles) equal to 2.")
	}
	y.mat <- y@data.mat
	y.vec <- y@data.vec
	R <- r*(r + 1)/2
	n <- sum(y.vec, na.rm = TRUE)
	r.i <- rowSums(y.mat, na.rm = TRUE)
	c.i <- colSums(y.mat, na.rm = TRUE)
	C <- choose(t + R - 1, R - 1)
	tables <- matrix(NA, nrow = R, ncol = C)
	indx <- 1
	for (i in 0:t) {
		for (j in 0:(t - i)) {
			tables[, indx] <- c(i, j, t - i - j)
			indx <- indx + 1
		}
	}
	const <- lfactorial(t) + lfactorial(t + R - 1) + lfactorial(2*n + r - 1) - lfactorial(2*t + r - 1) - lfactorial(n + t + R - 1) -
			(n - sum(diag(y.mat)))*log(2) - sum(lfactorial(r.i + c.i))
	bf.log <- vector(length = C)
	for (k in 1:C) {
		x.vec <- tables[, k]
		x <- matrix(NA, r, r)
		for (i in 1:r) {
			x[i, 1:i] <- x.vec[((i - 1)*i/2 + 1):(i*(i + 1)/2)]
		}
		z.i <- rowSums(x, na.rm = TRUE)
		d.i <- colSums(x, na.rm = TRUE)
		bf.log[k] <- const + (t - sum(diag(x)))*log(2) + sum(lfactorial(z.i + d.i)) - 2*sum(lfactorial(x.vec)) +
				sum(lfactorial(x.vec + y.vec))
	}
	out <- sum(exp(bf.log))
	return(out)
}

hwe.bf <- function(y) {
	if (class(y) != "HWEdata") {
		stop("y argument not of class 'HWEdata'. Type '?HWEdata' for help.")
	}
	y.mat <- y@data.mat
	r <- y@size
	R <- r*(r + 1)/2
	y.vec <- y@data.vec
	n <- sum(y.vec, na.rm = TRUE)
	r.i <- rowSums(y.mat, na.rm = TRUE)
	c.i <- colSums(y.mat, na.rm = TRUE)
	bf.log <- lfactorial(R - 1) + lfactorial(2*n + r - 1) - lfactorial(n + R - 1) - lfactorial(r - 1) - (n - sum(diag(y.mat)))*log(2) -
			sum(lfactorial(r.i + c.i)) + sum(lfactorial(y.vec))
	out <- exp(bf.log)
	return(out)
}

cip.tmp <- function(p11, p21, t, p) {
	R <- 3
	C <- choose(t + R - 1, R - 1)
	tables <- matrix(NA, nrow = R, ncol = C)
	indx <- 1
	for (i in 0:t) {
		for (j in 0:(t - i)) {
			tables[, indx] <- c(i, j, t - i - j)
			indx <- indx + 1
		}
	}
	const <- factorial(t)*factorial(t + 2)
	cip.terms <- numeric(C)
	pij <- c(p11, p21, 1 - p11 - p21)
	if (all(pij <= 1) & all(pij >= 0)) {
		for (k in 1:C) {
			x <- tables[, k]
			cip.terms[k] <- 2^x[2] * p^(2*x[1] + x[2]) * (1 - p)^(x[2] + 2*x[3]) * prod(pij^x) / prod(factorial(x))^2
		}
		out <- const*sum(cip.terms) }
	else { out <- NA }
	return(out)
}

cip.2 <- function(t, p, k = 30) {
	p11 <- p21 <- seq(0, 1, length.out = k)
	cip <- outer(X = p11, Y = p21, FUN = Vectorize(cip.tmp), t, p)
	cip[is.na(cip)] <- 0
	nrz <- nrow(cip)
	ncz <- ncol(cip)
	# Create a function interpolating colors in the range of specified colors
	colors <- colorRampPalette( c("gray","lightgray"))
	# Generate the desired number of colors from this palette
	nbcol <- 100
	color <- colors(nbcol)
	# Compute the z-value at the facet centres
	z <- cip
	zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
	# Recode facet z-values into color indices
	facetcol <- cut(zfacet, nbcol)
	persp(p11, p21, cip, xlab = "p_11", ylab = "p_21", zlab = "density", ticktype = "detailed", shade = .5, theta = -25, phi = 25,
			expand = 0.4, scale = TRUE, ltheta = 10, lphi = -10, cex.axis = 1.2, cex.lab = 1.2, col = color[facetcol])
}

ip.tmp <- function(p11, p21, t) {
	R <- 3
	C <- choose(t + R - 1, R - 1)
	tables <- matrix(NA, nrow = R, ncol = C)
	indx <- 1
	for (i in 0:t) {
		for (j in 0:(t - i)) {
			tables[, indx] <- c(i, j, t - i - j)
			indx <- indx + 1
		}
	}
	
	const <- factorial(t)*factorial(t + 2)/factorial(2*t + 1)
	ip.terms <- numeric(C)
	pij <- c(p11, p21, 1 - p11 - p21)
	if (all(pij <= 1) & all(pij >= 0)) {
		for (k in 1:C) {
			x <- tables[,k]
			ip.terms[k] <- 2^x[2] * factorial(2*x[1] + x[2]) * factorial(x[2] + 2*x[3]) * prod(pij^x) / prod(factorial(x))^2
		}
		out <- const*sum(ip.terms) }
	else { out <- NA }
	return(out)
}

ip.2 <- function(t, k = 30) {	
	p11 <- p21 <- seq(0, 1, length.out = k)
	ip <- outer(X = p11, Y = p21, FUN = Vectorize(ip.tmp), t)
	ip[is.na(ip)] <- 0
	nrz <- nrow(ip)
	ncz <- ncol(ip)
	# Create a function interpolating colors in the range of specified colors
	colors <- colorRampPalette( c("gray","lightgray"))
	# Generate the desired number of colors from this palette
	nbcol <- 100
	color <- colors(nbcol)
	# Compute the z-value at the facet centres
	z <- ip
	zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
	# Recode facet z-values into color indices
	facetcol <- cut(zfacet, nbcol)
	persp(p11, p21, ip, xlab = "p_11", ylab = "p_21", zlab = "density", ticktype = "detailed", shade = .5, theta = 50, phi = 35,
			expand = 0.4, scale = TRUE, ltheta = 40, lphi = -20, cex.axis = 1.2, cex.lab = 1.2, col = color[facetcol])
}

lik.multin <- function(y, p11, p21) {
	if (class(y) != "HWEdata") {
		stop("y argument not of class 'HWEdata'. Type '?HWEdata' for help.")
	}
	y.vec <- y@data.vec
	n <- sum(y.vec)
	out <- matrix(NA, nrow = length(p11), ncol = length(p21))
	for (i in 1:length(p11)) {
		for (j in 1:length(p21)) {
			prob <- c(p11[i], p21[j], 1 - p11[i] - p21[j])
			if (all(prob <= 1) & all(prob >= 0)) {
				out[i, j] <- dmultinom(y.vec, prob = prob) }
			else{ out[i, j] <- NA }
		}
	}
	const <- lfactorial(n) + sum(lgamma(y.vec + 1)) - sum(lfactorial(y.vec)) - lgamma(n + length(y.vec))
	const <- exp(const)
	out <- out/const
	return(out)
}

hwe.ibf.plot <- function(y, t.vec, M = 100000, bf = FALSE) {
	if (class(y) != "HWEdata") {
		stop("y argument not of class 'HWEdata'. Type '?HWEdata' for help.")
	}
	n <- sum(y@data.vec, na.rm = TRUE)
	H <- length(t.vec)
	out.mc <- matrix(NA, nrow = H, ncol = 2)
	for (h in 1:H) {
		print(paste(h, "/", H, sep = ""), quote = FALSE)
		res <- hwe.ibf.mc(y, t.vec[h], M, verbose = FALSE)
		out.mc[h,] <- c(log10(res@bf), res@npp)
	}
	if (y@size == 2) {
		out.exact <- matrix(NA, nrow = H, ncol = 2)
		hwe.ibf.exact <- Vectorize(hwe.ibf, "t")
		res <- hwe.ibf.exact(t.vec, y = y)
		out.exact <- cbind(log10(res), 1/(1 + res))
	}
	res <- hwe.bf(y)
	out.std <- c(log10(res), 1/(1 + res))
	if (bf) {
		idx = 1
		ytxt = "Bayes Factor of Alternative vs. Null (log10)"
	}
	else {
		idx = 2
		ytxt = "Null Posterior Probability"
	}
	if (y@size != 2) {
		out <- list(mc = out.mc, std = out.std)
	}
	else {
		out <- list(mc = out.mc, std = out.std, exact = out.exact)
	}
	plot(t.vec/n, out.mc[, idx], type = "n", ylab = ytxt, xlab = "f = t/n", main = "", cex.axis = 1.75,
			cex.lab = 1.75, col = gray(.5), ylim = c(min(out.std[idx], out.mc[, idx]), max(out.std[idx], out.mc[, idx])))
	abline(h = out.std[idx], col = gray(.5), lty = "longdash")
	if (y@size == 2) lines(t.vec/n, out.exact[, idx], col = gray(.3), lty = "longdash", lwd = 2)
	lines(smooth.spline(t.vec/n, out.mc[, idx]), col = gray(0), lwd = 2)
	invisible(out)
}
