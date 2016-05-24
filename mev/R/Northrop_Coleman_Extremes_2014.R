#######################################################################################################
# Functions to create plots and tests in article
# Northrop, Paul J and Claire L Coleman (2014), Extremes
# "Improved threshold diagnostic plots for extreme value analyses"
#
#######################################################################################################
#------------------------------------------------------------------------------#
#                             Main function                                    #
#------------------------------------------------------------------------------#
#' Score and likelihood ratio tests fit of equality of shape over multiple thresholds
#'
#' The function returns a P-value path for the score testand/or likelihood ratio
#' test for equality of the shape parameters over
#' multiple thresholds under the generalized Pareto model.
#'
#' @param x  raw data
#' @param u \code{m}-vector of thresholds (sorted from smallest to largest)
#' @param GP.fit function used to optimize the generalized Pareto model.
#' @param do.LRT boolean indicating whether to perform the likelihood ratio test (in addition to the score test)
#' @param size level at which a horizontal line is drawn on multiple threshold plot
#' @param my.xlab (optional) x-axis label
#' @param xi.tol numerical tolerance for threshold distance; if the absolute value of \code{xi1.hat} is less than \code{xi.tol} use linear interpolation
#'                to evaluate score vectors, expected Fisher information matrices, Hessians
#'
#' @details The default method is \code{"Grimshaw"} using the reduction of the parameters to a one-dimensional
#' maximization. Other options are one-dimensional maximization of the profile the \code{nlm} function or \code{optim}.
#' Two-dimensional optimisation using 2D-optimization \code{\link[ismev]{ismev}} using the routine
#' from \code{gpd.fit} from the \code{ismev} library, with the addition of the algebraic gradient.
#' The choice of \code{GP.fit} should make no difference but the options were kept.
#' \bold{Warning}: the function is not robust
#' and will not recover from failure of the maximization routine, returning various error messages.
#'
#'
#' @references Grimshaw (1993). Computing Maximum Likelihood Estimates for the Generalized
#'  Pareto Distribution, \emph{Technometrics}, \bold{35}(2), 185--191.
#' @references Northrop & Coleman (2014). Improved threshold diagnostic plots for extreme value
#' analyses, \emph{Extremes}, \bold{17}(2), 289--303.
#' @references Wadsworth & Tawn (2012). Likelihood-based procedures for threshold
#' diagnostics and uncertainty in extreme value modelling, \emph{J. R. Statist. Soc. B}, \bold{74}(3), 543-???567.
#' @export
#' @author Paul J. Northrop and Claire L. Coleman
#' @return a plot of P-values for the test at the different thresholds \code{u}
#' @examples
#' \dontrun{
#' library(ismev)
#' data(rain)
#' u <- quantile(rain, seq(0.85,0.99,by=0.01))
#' NC.diag(rain, u, size=0.05)
#' }
NC.diag <- function(x, u, GP.fit = c("Grimshaw", "nlm", "optim", "ismev"),
			do.LRT = FALSE, size = NULL, my.xlab = NULL, xi.tol = 1e-3) {
		if (any(diff(u) <= 0)) {
			warning("Thresholds supplied in u are not in increasing order")
		}
		u <- sort(u)
		n_u <- length(u)      # total number of thresholds
		#------------------------------------------------------------------------------#
		# 1. Fit GP distribution to excesses of u[i], i=1, ..., n_u                      #
		#------------------------------------------------------------------------------#
		GP.fit <-	match.arg(arg = GP.fit, choices = c("Grimshaw", "nlm", "optim", "ismev"))
		z <- list()                 # list to store the results
		z$thresh <- u               # all thresholds
		z$nexc <-	unlist(lapply(u, function(y) {
				sum(x > y)
			})) # number of excesses of each threshold
		z$n.between <-	c(-diff(z$nexc), z$nexc[n_u])
		# sample sizes between thresholds (and above the highest threshold)
		for (i in 1:n_u) {
			# loop over all thresholds
			if (GP.fit == "nlm") {
				temp <-	.gpd_1D_fit(x, u[i], show = F, xi.tol = xi.tol, calc.se = F) # threshold u[j1], 1D max, algebraic Hessian
				phi.init <-	temp$mle[2] / temp$mle[1] # better initial estimate of phi
				temp <-	.GP_1D_fit_nlm(x, u[i], init.val = phi.init, gradtol = 1e-20, steptol = 1e-20, calc.se = F)
				# 1D max, use nlm to get gradients v close to zero
			} else if (GP.fit == "ismev") {
				temp <-	.gpd_2D_fit(x, u[i], show = F) # threshold u[i], ismev GP fitting function
				temp <-	.gpd_2D_fit(x, u[i], show = F, siginit = temp$mle[1], shinit = temp$mle[2], method = "BFGS", reltol =	1e-30, abstol = 1e-30)
			} else if (GP.fit == "optim") {
				temp <-	.gpd_1D_fit(x, u[i], show = F, xi.tol = xi.tol)					# threshold u[i], 1D max, algebraic Hessian
				temp <-	.gpd_1D_fit(x, u[i], show = F, xi.tol = xi.tol, phi.input = temp$mle[2] / temp$mle[1], reltol = 1e-30, abstol = 1e-30) # threshold u[i], 1D max, algebraic Hessian
			} else if (GP.fit == "Grimshaw") {
				yy <- x[x > u[i]] - u[i] # thresholds excesses
				pjn <-	.gpd_grimshaw(yy)  # Grimshaw (1993) function, note: k is -xi, a is sigma
				temp <- list()
				temp$mle <- c(pjn$a, -pjn$k)        # mle for (sigma, xi)
				sc <- rep(temp$mle[1], length(yy));
				xi <- temp$mle[2]
				temp$nllh <- sum(log(sc)) + sum(log(1 + xi * yy / sc) * (1 / xi + 1))
			}
			z$xi.mle[i]    <- temp$mle[2]    # MLE of xi
			z$sigma.mle[i] <- temp$mle[1]    # MLE of sigma1
			z$nllh[i] <- temp$nllh           # negated log-likelihood at MLE
		}
		#.....................# end of loop over thresholds
		#------------------------------------------------------------------------------#
		# 2. For each threshold u[i], i=1, ..., m-1 calculate score-based p-value for  #
		#    the test of H_0: xi[i]=...=xi[m], equivalently H_0: phi_1=...=phi_m       #
		#    ... and produce plot of the p-values vs threshold.                        #
		#------------------------------------------------------------------------------#
		# Create functions score.test(), mult.u.gpd.fit() and mult.thresh.LR.test() ...
		##########################################################################################
		#-------------------------- Start of function score.test() ------------------------------#
		##########################################################################################
		score.test <- function(my.data, m, vi, wi, pars, sigma1, xi1) {
			#--------------------- if xi.hat is not very close to 0 ... ---------------------------#
			if (abs(xi1) >= xi.tol) {
				score <-	.score_algebraic(my.data, pars, wi, vi, m)    # score vector under H_0
				e.info <-	n.exc * .exp_info_algebraic(pars, wi, vi, m)   # expected information under H_0
			}
			#----------------------- if xi.hat is very close to 0 ... -----------------------------#
			if (abs(xi1) < xi.tol) {
				delta <- 2 * xi.tol   # evaluate score/info at xi+delta and xi-delta
				phis <- 1 / (sigma1 / (xi1 + delta) + cumsum(c(0, wi[-m])))     # values of phi under H_0
				pars1 <- c(sigma1, phis)  # estimates of sigma1, phi_1, ..., phi_m under H_0.
				phis <-	1 / (sigma1 / (xi1 - delta) + cumsum(c(0, wi[-m])))     # values of phi under H_0
				pars2 <-	c(sigma1, phis)  # estimates of sigma1, phi_1, ..., phi_m under H_0.
				score1 <-	.score_algebraic(my.data, pars1, wi, vi, m)       # score vector at xi+delta
				score2 <-	.score_algebraic(my.data, pars2, wi, vi, m)       # score vector at xi-delta
				score <-	(score1 + score2) / 2     # average of the two scores
				info1 <- n.exc * .exp_info_algebraic(pars1, wi, vi, m) # expected information at xi+delta
				info2 <- n.exc * .exp_info_algebraic(pars2, wi, vi, m)  # expected information at xi-delta
				e.info <-	(info1 + info2) / 2                        # average of the two information matrices
			}
			#.......... Start of test.stat.calc() ...........
			test.stat.calc <- function(score, info) {
				my.svd <- 	svd(info) # SVD of information matrix
				vec <- t(score) %*% my.svd$v  # avoid matrix inversion
				stat <-	vec %*% diag(1 / my.svd$d, nrow = m + 1) %*% t(vec)   # score test statistic
				stat
			}
			#........... End of test.stat.calc() ............

			e.stat <- 	test.stat.calc(score, e.info)   # score test statistic
			e.p <- pchisq(e.stat, df = m - 1, lower.tail = F)   # p-value
			c(e.stat, e.p)
		}
		#.......................................................
		# end of function score.test

		##########################################################################################
		#-------------------------- End of function score.test() --------------------------------#
		##########################################################################################

		##########################################################################################
		#----------------------- Start of function mult.u.gpd.fit() -----------------------------#
		##########################################################################################

		mult.u.gpd.fit <-	function(y, m, v, w, npy = 365, method = "Nelder-Mead", maxit = 10000,
			init.ests =	NULL, ...) {
				negated.mult.log.likelihood <- function(sig) {
					# sig: (sigma1, xi_1, ..., xi_m)
					# y: excesses of threshold
					sigma1 <- sig[1] # sigma_1
					if (sigma1 <= 0)
						return(1e30) # Need sigma_1 > 0
					xi <-	sig[-1]   # (xi_1, ..., xi_m)
					sigma <- sigma1 + cumsum(c(0, xi[-m] * w[-m]))  # (sigma_1, ..., sigma_m)
					phi <-	xi / sigma     # (phi_1, ..., phi_m)
					if (any(1 + phi[-m] * w[-m] <= 0)){
						return(1e30)        # Need all elements of 1+phi*w/sigma > 0
					}
					Ij <-	unlist(lapply(y, function(sig){
							sum(sig - v > 0)})) # interval indicators
					if (any(1 + phi[Ij] * (y - v[Ij]) <= 0))
						return(1e30)    # Need all elements of 1+phi[Ij]*(y-v[Ij]) > 0
					aj <-	c(0, cumsum(log(1 + phi[-m] * w[-m]) / sigma[-m] / phi[-m])) # -log(p_j), j=1, ..., m
					pj <-	exp(-aj)    # P(Y > v_j), j=1, ..., m
					bj <- log(sigma)
					dj <- log(1 + phi[Ij] * (y - v[Ij]))
					ej <- log(1 + phi[Ij] * (y - v[Ij])) / sigma[Ij] / phi[Ij]
					sum(aj[Ij] + bj[Ij] + dj + ej)
				}

				fscale <- negated.mult.log.likelihood(init.ests)
				temp <- optim(init.ests, negated.mult.log.likelihood, hessian = F, method = method,
						control = list(maxit = maxit, fnscale = fscale, ...)
					)
				zz <- list()
				zz$mle <- temp$par
				zz$nllh <- temp$value
				zz$conv <- temp$convergence
				zz$counts <- temp$counts
				zz$message <- temp$message
				invisible(zz)
			}

		##########################################################################################
		#------------------------- End of function mult.u.gpd.fit() -----------------------------#
		##########################################################################################

		##########################################################################################
		#--------------------- Start of function mult.thresh.LR.test() --------------------------#
		##########################################################################################

		mult.thresh.LR.test <-function(my.data, m, vi, wi, init.ests = NULL, null.nllh = NULL) {
				nllh2 <- mult.u.gpd.fit(my.data, m, vi, wi, init.ests = init.ests)$nllh
				LRT.stat <-	ifelse(is.null(null.nllh), 2 * (z$nllh[i] - nllh2), 2 * (null.nllh - nllh2))
				pvalue <- pchisq(LRT.stat, m - 1, lower.tail = FALSE)
				c(LRT.stat, pvalue)
			}

		##########################################################################################
		#----------------------- End of function mult.thresh.LR.test() --------------------------#
		##########################################################################################


		for (i in 1:(n_u - 1)) {
			# loop from u[1] to u[n_u-1]
			excess.data <- x[x > u[i]] - u[i] # excesses of threshold u[i]
			n.exc <-	z$nexc[i]   # number of excesses of current threshold
			# (calculated at top of function)
			m <- n_u - i + 1 # number of thresholds at or above u[i]
			z$df[i] <-	m - 1 # df of null chi-squared distribution
			vi <-	u[i:n_u] - u[i] # thresholds u[i], ..., u[n_u] relative to u[i]
			wi <-	c(diff(vi), NA) # differences between thresholds (wi[m] not used)
			phis <-	1 / (z$sigma.mle[i] / z$xi.mle[i] + cumsum(c(0, wi[-m]))) # values of phi under H_0
			pars <-	c(z$sigma.mle[i], phis) # estimates of sigma1, phi_1, ..., phi_m under H_0.
			# MLEs of sigma1 and xi1 under H_0, for passing to score.test() ...
			sigma1 <- z$sigma.mle[i];  xi1 <- z$xi.mle[i]
			temp <-	score.test(excess.data, m, vi, wi, pars, sigma1, xi1) # do score test
			z$e.test.stats[i] <- temp[1]
			z$e.p.values[i] <- temp[2]
			if (do.LRT) {
				null.init.ests <- c(sigma1, rep(xi1, m))
				temp <-	mult.thresh.LR.test(excess.data, m, vi, wi, init.ests = null.init.ests) # do LR test
				z$LRT.p.values[i] <- temp[2]
				z$LRT.test.stats[i] <- temp[1]
			}
		}#.......................# end of loop over thresholds
		z$u <-	u[1:(n_u - 1)]   # (lowest) thresholds for each test

		# Produce the plot ......

		my.ylab <- "p-value";
		if (is.null(my.xlab))
			my.xlab <- "threshold"
		plot(z$u, z$e.p.values, type = "b", xlab = my.xlab, ylab = "p-value",
			pch = 16, ylim = c(0, 1)
		)
		#  axis(3, at=u[1:num.u], labels=z$n.between, cex.axis=0.7)
		axis(3, at = u[1:n_u], labels = z$nexc[1:n_u], cex.axis = 0.7)
		if (do.LRT)
			lines(z$u, z$LRT.p.values, type = "b", lty = 4, pch = 2)
		if (!is.null(size))
			abline(h = size, lty = 2)
		#
		invisible(z)
	}

#------------------------------------------------------------------------------#
#                  Algebraic calculation of score vector                       #
#------------------------------------------------------------------------------#

#' Algebraic score
#'
#' @param y vector of excesses of lowest threshold \code{u1}
#' @param x parameter vector (\code{sigma1}, \code{phi_1}, \ldots, \code{phi_m})
#' @param v thresholds relative to lowest threshold
#' @param w differences between thresholds (\code{w[m]} not used)
#' @param m number of thresholds
.score_algebraic <- function(y, x, w, v, m) {
	n <- length(y)                                # sample size
	sigma1 <- x[1]                                # sigma_1
	phi <- x[-1]                                  # (phi_1, ..., phi_m)
	sigma <- sigma1 * c(1, cumprod(1 + phi[-m] * w[-m])) # (sigma_1, ..., sigma_m)
	xi <- sigma * phi                               # (xi_1, ..., xi_m)
	Ij <- unlist(lapply(y, function(x){
		sum((x - v) > 0)})) # interval indicators
	h <- phi * sigma / sigma1                         # (h_1, ..., h_m)
	##################### Derivatives of bj ...
	db.s1 <- n / sigma1
	db.phi <-	c(unlist(lapply(1:(m - 1), function(x){
			sum(Ij > x) * w[x] / (1 + phi[x] * w[x])})), 0)
	##################### Derivatives of dj ...
	dd.s1 <- 0
	dd.phi <-	unlist(lapply(1:m, function(x){
			sum((y[Ij == x] - v[x]) / (1 + phi[x] * (y[Ij == x] - v[x])))
		}))
	##################### Derivatives of ej ...
	de.s1 <- sum(-sigma1^(-2) * h[Ij]^(-1) * log(1 + phi[Ij] * (y - v[Ij])))
	de.phi.fn <- function(x) {
	temp1 <-	h[x]^(-1) * (-phi[x]^(-1) * log(1 + phi[x] * (y[Ij == x] - v[x])) +
					(y[Ij == x] - v[x]) * (1 + phi[x] * (y[Ij == x] - v[x]))^(-1))
		if (x == m){
			return(sigma1^(-1) * sum(temp1))
		}
		ind1 <- Ij %in% (x + 1):m
		ind2 <- Ij[ind1]
		temp2 <- -w[x] * (1 + phi[x] * w[x])^(-1) * (h[Ij][ind1]^(-1) * log(1 + phi[Ij][ind1] *
					(y[ind1] - v[Ij][ind1])))
		sigma1^(-1) * (sum(temp1) + sum(temp2))
	}
	de.phi <- unlist(lapply(1:m, de.phi.fn))
	##################### Derivatives of aj ...
	aj <- 	c(0, cumsum(log(1 + phi[-m] * w[-m]) / sigma[-m] / phi[-m])) # -log(p_j), j=1, ..., m
	da.s1 <- -sigma1^(-1) * sum(aj[Ij])
	da.phi.fn <- function(x) {
		n.Ijs <- table(c(Ij, 1:m)) - 1 # frequencies of values of Ij
		n.Ijs <- n.Ijs[-1]
		temp1 <- h[x]^(-1) * (-phi[x]^(-1) * log(1 + phi[x] * w[x]) + w[x] * (1 + phi[x] *	w[x])^(-1)) # B(x, x)
		my.zeros <- numeric(m - 1)
		my.zeros[x] <- temp1
		temp1 <- my.zeros
		n.Ijs <- as.numeric(n.Ijs)
		if (x == (m - 1)){
			return(sigma1^(-1) * sum(temp1 * n.Ijs))
		}
		if (x < (m - 1)) {
			ind <- (x + 1):(m - 1)
			temp2 <- -w[x] * (1 + phi[x] * w[x])^(-1) * h[ind]^(-1) * log(1 + phi[ind] *	w[ind])
		}
		temp1[(x + 1):(m - 1)] <- temp2
		sigma1^(-1) * sum(cumsum(temp1) * n.Ijs)
	}
	#
	da.phi <- c(unlist(lapply(1:(m - 1), da.phi.fn)), 0)
	########### Return score vector ...
	-c(da.s1 + db.s1 + dd.s1 + de.s1, da.phi + db.phi + dd.phi + de.phi)
}

#------------------------------------------------------------------------------#
#                  Algebraic calculation of expected information               #
#------------------------------------------------------------------------------#

#' Algebraic calculation of the expected information
#'
#' @param x parameter vector: (\code{sigma1}, \code{phi_1}, \ldots, \code{phi_m})
#' @param v thresholds relative to lowest threshold
#' @param v thresholds relative to lowest threshold
#' @param w differences between thresholds (\code{w[m]} not used)
#' @param m number of thresholds
#'
.exp_info_algebraic <- function(x, w, v, m) {
	sigma1 <- x[1]                                  # sigma_1
	phi <- x[-1]                                    # (phi_1, ..., phi_m)
	sigma <-	sigma1 * c(1, cumprod(1 + phi[-m] * w[-m]))   # (sigma_1, ..., sigma_m)
	xi <- sigma * phi     # (xi_1, ..., xi_m)
	aj <-	c(0, cumsum(log(1 + phi[-m] * w[-m]) / sigma[-m] / phi[-m])) # -log(p_j), j=1, ..., m
	pj <-	exp(-aj)     # P(Y > v_j), j=1, ..., m
	qj <-	c(pj[-m] * (1 - (1 + phi[-m] * w[-m])^(-1 / xi[-m])), pj[m]) # P(v_j < Y < v_{j+1}), j=1, ..., m
	h <- phi * sigma / sigma1   # (h_1, ..., h_m)
	#
	##################### Various integrals ...
	#
	I0b <- function(b, j) {
		t1 <- (1 + b * xi[j])^(-1)
		if (j == m)
			return(t1)
		t1 - (1 + phi[j] * w[j])^(-b - 1 / xi[j]) * t1
	}
	#
	I1b <- function(b, j) {
		t1 <- sigma[j] / (1 + (b - 1) * xi[j]) / (1 + b * xi[j])
		if (j == m)
			return(t1)
		t2 <- (1 + phi[j] * w[j])^(-b - 1 / xi[j]) / (1 + b * xi[j])
		t3 <- (1 + phi[j] * w[j])^(1 - b - 1 / xi[j]) / (1 + (b - 1) * xi[j])
		t1 + (t2 - t3) / phi[j]
	}
	#
	I2b <- function(b, j) {
		t1 <- 2 * sigma[j]^2 / (1 + (b - 2) * xi[j]) / (1 + (b - 1) * xi[j]) /
			(1 + b * xi[j])
		if (j == m)
			return(t1)
		t2 <- (1 + phi[j] * w[j])^(2 - b - 1 / xi[j]) / (1 + (b - 2) * xi[j])
		t3 <- (1 + phi[j] * w[j])^(1 - b - 1 / xi[j]) / (1 + (b - 1) * xi[j])
		t4 <- (1 + phi[j] * w[j])^(-b - 1 / xi[j]) / (1 + b * xi[j])
		t1 - (t2 - 2 * t3 + t4) / phi[j]^2
	}
	#
	J <- function(j) {
		t1 <- xi[j]
		if (j == m)
			return(t1)
		t2 <- -(1 + phi[j] * w[j])^(-1 / xi[j]) * log(1 + phi[j] * w[j])
		t3 <- -xi[j] * (1 + phi[j] * w[j])^(-1 / xi[j])
		t1 + t2 + t3
	}
	#
	##################### Derivatives of bj (constant w.r.t. y) ...
	#
	# Note: all contributions are zero apart from the diagonal elements
	#       for sigma1, phi_1, ..., phi_{m-1}
	b.exp.info <-	matrix(0, m + 1, m + 1)              # matrix for expected information from b
	db.phi.phi <-	c(unlist(lapply(1:(m - 1), function(x){
			- w[x]^2 * (1 + phi[x] * w[x])^(-2) * sum(qj[(x + 1):m])})), 0)
	diag(b.exp.info) <- c(-1 / sigma1^2, db.phi.phi)
	#
	##################### Derivatives of dj ...
	# Note: all contributions are zero apart from the diagonal elements
	#       for phi_1, ..., phi_m.
	d.exp.info <- matrix(0, m + 1, m + 1)              # matrix for expected information from b
	dd.phi.phi <- unlist(lapply(1:m, function(x){
		- pj[x] * I2b(2, x)}))
	diag(d.exp.info) <- c(0, dd.phi.phi)
	##################### Derivatives of ej ...
	#
	e.exp.info <-	matrix(0, m + 1, m + 1)              # matrix for expected information from b
	de.s1.s1 <- 2 * sigma1^(-3) * sum(h^(-1) * pj * unlist(lapply(1:m, J)))
	#
	de.phi.fn <- function(x) {
		temp1 <- h[x]^(-1) * (-phi[x]^(-1) * J(x) + I1b(1, x))
		if (x == m)
			return(sigma1^(-1) * temp1 * pj[x])
		temp2 <-	-w[x] * (1 + phi[x] * w[x])^(-1) * h[(x + 1):m]^(-1) * unlist(lapply((x +	1):m, J))
		sigma1^(-1) * (temp1 * pj[x] + sum(temp2 * pj[(x + 1):m]))
	}
	#
	de.phi <- unlist(lapply(1:m, de.phi.fn))
	#
	de.phi.phi.fn <- function(x) {
		temp1 <-	h[x]^(-1) * (2 * phi[x]^(-2) * J(x) - 2 * phi[x]^(-1) * I1b(1, x) -	I2b(2, x))
		if (x == m)
			return(sigma1^(-1) * temp1 * pj[x])
		temp2 <-
			2 * w[x]^2 * (1 + phi[x] * w[x])^(-2) * h[(x + 1):m]^(-1) * unlist(lapply((x +	1):m, J))
		sigma1^(-1) * (temp1 * pj[x] + sum(temp2 * pj[(x + 1):m]))
	}
	#
	de.phi.phi <- unlist(lapply(1:m, de.phi.phi.fn))
	#
	de.phil.phik.fn <- function(x) {
		k <- x[1]; l <- x[2]
		temp1 <-		h[k]^(-1) * w[l] * (1 + w[l] * phi[l])^(-1) * (phi[k]^(-1) * J(k) -	I1b(1, k))
		if (k == m)
			return(sigma1^(-1) * temp1 * pj[k])
		temp2 <-
			w[k] * (1 + phi[k] * w[k])^(-1) * w[l] * (1 + phi[l] * w[l])^(-1) *
			h[(k + 1):m]^(-1) * unlist(lapply((k + 1):m, J))
		sigma1^(-1) * (temp1 * pj[k] + sum(temp2 * pj[(k + 1):m]))
	}
	k.vals <- rep(2:m, times = 2:m - 1) # k > l
	l.vals <- unlist(lapply(2:m - 1, function(x)
		seq(from = 1, to = x)))
	kl.vals <- cbind(k.vals, l.vals)
	rev.kl.vals <- kl.vals[, 2:1]
	if (m == 2) {
		kl.vals <- matrix(kl.vals, nrow = 1, ncol = 2) # if m=2 make kl.vals a matrix
		rev.kl.vals <- matrix(kl.vals[, 2:1], nrow = 1, ncol = 2) # if m=2 make rev.kl.vals a matrix
	}
	de.phil.phik <- unlist(apply(kl.vals, 1, de.phil.phik.fn))
	#
	diag(e.exp.info) <- c(de.s1.s1, de.phi.phi)
	e.exp.info[1, 2:(m + 1)] <-	e.exp.info[2:(m + 1), 1] <- -sigma1^(-1) * de.phi
	e.exp.info[kl.vals + 1] <-	e.exp.info[rev.kl.vals + 1] <- de.phil.phik
	##################### Derivatives of aj ...
	a.exp.info <- 	matrix(0, m + 1, m + 1)              # matrix for expected information from b
	aj <-	c(0, cumsum(log(1 + phi[-m] * w[-m]) / sigma[-m] / phi[-m])) # -log(p_j), j=1, ..., m
	da.s1.s1 <- 2 * sigma1^(-2) * sum(aj * qj)
		da.phi.fn <- function(x) {
		temp1 <-	h[x]^(-1) * (-phi[x]^(-1) * log(1 + phi[x] * w[x]) + w[x] * (1 + phi[x] *	w[x])^(-1)) # B(x, x)
		my.zeros <- numeric(m - 1)
		my.zeros[x] <- temp1
		temp1 <- my.zeros
		if (x == (m - 1))
			return(sigma1^(-1) * sum(temp1 * qj[-1]))
		if (x < (m - 1)) {
			ind <- (x + 1):(m - 1)
			temp2 <- -w[x] * (1 + phi[x] * w[x])^(-1) * h[ind]^(-1) * log(1 + phi[ind] *	w[ind])
		}
		temp1[(x + 1):(m - 1)] <- temp2
		sigma1^(-1) * sum(cumsum(temp1) * qj[-1])
	}
	#
	da.phi <- c(unlist(lapply(1:(m - 1), da.phi.fn)), 0)
	da.s1.phi <- -sigma1^(-1) * da.phi
	#
	da.phi.phi.fn <- function(x) {
		y.1v <- 1 + phi[x] * w[x]
		y.v <- w[x]
		temp1 <-	h[x]^(-1) * (2 * phi[x]^(-2) * log(y.1v) - 2 * phi[x]^(-1) * y.v *
					y.1v^(-1) - y.v^2 * y.1v^(-2))
		my.zeros <- numeric(m - 1)
		my.zeros[x] <- temp1
		temp1 <- my.zeros
		if (x == (m - 1))
			return(sigma1^(-1) * sum(temp1 * qj[-1]))
		if (x < (m - 1)) {
			ind <- (x + 1):(m - 1)
			temp2 <-	2 * w[x]^2 * (1 + phi[x] * w[x])^(-2) * h[ind]^(-1) *
				log(1 + phi[ind] *	w[ind])
		}
		temp1[(x + 1):(m - 1)] <- temp2
		sigma1^(-1) * sum(cumsum(temp1) * qj[-1])
	}
	#
	da.phi.phi <- c(unlist(lapply(1:(m - 1), da.phi.phi.fn)), 0)
	#
	da.phil.phik.fn <- function(x) {
		k <- x[1]; l <- x[2]
		y.1v <- 1 + phi[k] * w[k]
		y.v <- w[k]
		temp1 <-	h[k]^(-1) * w[l] * (1 + w[l] * phi[l])^(-1) * (phi[k]^(-1) * log(y.1v) -
					y.v * y.1v^(-1))
		my.zeros <- numeric(m - 1)
		my.zeros[k] <- temp1
		temp1 <- my.zeros
		if (k == (m - 1))
			return(sigma1^(-1) * sum(temp1 * qj[-1]))
		if (k < (m - 1)) {
			ind <- (k + 1):(m - 1)
			temp2 <-	w[k] * (1 + phi[k] * w[k])^(-1) * w[l] * (1 + phi[l] * w[l])^(-1) *
				h[ind]^(-1) * log(1 + phi[ind] * w[ind])
		}
		temp1[(k + 1):(m - 1)] <- temp2
		sigma1^(-1) * sum(cumsum(temp1) * qj[-1])
	}
	diag(a.exp.info) <- c(da.s1.s1, da.phi.phi)
	a.exp.info[1, 2:(m + 1)] <- a.exp.info[2:(m + 1), 1] <- da.s1.phi
	if (m > 2) {
		k.vals <-	rep(2:(m - 1), times = 1:(m - 2))  # k > l
		l.vals <- unlist(lapply(1:(m - 2), function(x){
			seq(from = 1, to = x)}))
		kl.vals <- cbind(k.vals, l.vals)
		rev.kl.vals <- kl.vals[, 2:1]
		if (m == 3) {
			kl.vals <-	matrix(kl.vals, nrow = 1, ncol = 2) # if m=3 make kl.vals a matrix
			rev.kl.vals <- matrix(kl.vals[, 2:1], nrow = 1, ncol = 2) # if m=3 make rev.kl.vals a matrix
		}
		da.phil.phik <- unlist(apply(kl.vals, 1, da.phil.phik.fn))
		a.exp.info[kl.vals + 1] <-	a.exp.info[rev.kl.vals + 1] <- da.phil.phik
	}
	#
	##################### Return observed information ...
	exp.info <- a.exp.info + b.exp.info + d.exp.info + e.exp.info
	#
	exp.info
}


#------------------------------------------------------------------------------#
#               Fit GP (sigma, xi) distribution using 1D optimisation            #
#------------------------------------------------------------------------------#
.gpd_1D_fit <-	function (xdat, threshold, npy = 365, ydat = NULL, sigl = NULL,
		shl = NULL, siglink = identity, shlink = identity, siginit = NULL,
		shinit = NULL, show = TRUE, method = "BFGS", maxit = 10000, xi.tol = 1e-3,
	  phi.input =	NULL, 	calc.se = F, ...) {
		z <- list()
		npsc <- length(sigl) + 1
		npsh <- length(shl) + 1
		n <- length(xdat)
		z$trans <- FALSE
		if (is.function(threshold))
			stop("`threshold' cannot be a function")
		u <- rep(threshold, length.out = n)
		if (length(unique(u)) > 1)
			z$trans <- TRUE
		xdatu <- xdat[xdat > u]
		xind <- (1:n)[xdat > u]
		u <- u[xind]
		ex.data <- xdatu - u
		in2 <- sqrt(6 * var(xdat)) / pi
		in1 <- mean(xdat, na.rm = TRUE) - 0.57722 * in2
		phi.init <- (2^0.1 - 1) / median(xdatu - u)

		if (is.null(sigl)) {
			sigmat <- as.matrix(rep(1, length(xdatu)))
			if (is.null(siginit)){
				siginit <- in2
		}
		} else {
			z$trans <- TRUE
			sigmat <- cbind(rep(1, length(xdatu)), ydat[xind, sigl])
			if (is.null(siginit))
				siginit <- c(in2, rep(0, length(sigl)))
		}
		if (is.null(shl)) {
			shmat <- as.matrix(rep(1, length(xdatu)))
			if (is.null(shinit)) 	shinit <- 0.1
		} else {
			z$trans <- TRUE
			shmat <- cbind(rep(1, length(xdatu)), ydat[xind, shl])
			if (is.null(shinit))
				shinit <- c(0.1, rep(0, length(shl)))
		}
		init <- c(siginit, shinit)
		z$model <- list(sigl, shl)
		z$link <- deparse(substitute(c(siglink, shlink)))
		z$threshold <- threshold
		z$nexc <- length(xdatu)
		z$data <- xdatu

		GP.1D.negloglik <- function(phi) {
			# negated 1D loglikelihood
			k <- length(ex.data)
			if (phi == 0) {
				sigma.mle <- (1 / k) * sum(ex.data)
				return(k * log(sigma.mle) + (1 / sigma.mle) * sum(ex.data))
			}
			zz <- 1 + phi * ex.data
			if (min(zz) <= 0)
				return(1e30)
			xi.of.phi <- mean(log(zz))
			sigma.of.phi <- xi.of.phi / phi
			if (sigma.of.phi <= 0)
				return(1e30)
			k * log(sigma.of.phi) + (1 + 1 / xi.of.phi) * sum(log(1 + phi * ex.data))
		}#................................# end of GP.1D.negloglik()
		#
		gp.1D.grad <- function(a) {
			# gradient of negated log-likelihood
			phi <- a
			n <- length(ex.data)
			yy <- 1 + phi * ex.data
			xi.phi <- mean(log(yy))
			- n / phi + (1 + 1 / xi.phi) * sum(ex.data / yy)
		}
		#
		if (!is.null(phi.input)){	phi.init <- phi.input}
		temp <- 	try(optim(phi.init, GP.1D.negloglik, gr = gp.1D.grad, hessian = FALSE, method = method,
				control = list(maxit = maxit, ...)
			))
		if(is.character(temp)){
		 z$conv <- 50
		 return(z);
		}
		phi <- temp$par
		zz <- 1 + phi * (xdatu - u)
		if (min(zz) <= 0){
				z$conv <- 50
				return(z)
		}
		xi <- mean(log(zz))
		sc <- xi / phi
		z$mle <- c(sc, xi)
		#
		if (calc.se) {
			if (abs(xi) >= xi.tol && xi > -0.5)
				z$cov <- solve(.gpd_obs_info(scale=sc, shape=xi, data=xdatu - u))
			if (abs(xi) < xi.tol  && xi > -0.5) {
				delta <- 2 * xi.tol # evaluate observed information at xi+delta and xi-delta
				o.info1 <- .gpd_obs_info(scale=sc, shape=xi + delta, data=xdatu - u)
				o.info2 <- .gpd_obs_info(scale=sc, shape=xi - delta, data=xdatu - u)
				z$cov <- solve((o.info1 + o.info2) / 2)
			}
		}
		#
		z$conv <- temp$convergence
		z$counts <- temp$counts
		z$nllh <- temp$value
		z$vals <- cbind(sc, xi, u)
		z$rate <- length(xdatu) / n
		if (calc.se){
		#	z$se <- tryCatch(sqrt(diag(z$cov)), error = function(e) NULL)
		}
		z$n <- n
		z$npy <- npy
		z$xdata <- xdat
		if (show)
			print(z[c(4, 5, 9, 10, 7, 12, 13)])
				uu <- unique(u)
		z$mle.t <-	c(z$mle[1] - z$mle[2] * uu, z$mle[2]) # sigma-xi*u replaces sigma

		if (calc.se) {
			d <-	matrix(c(1, -uu, 0, 1), 2, 2, byrow = T)       # derivatives of sigma-xi*u and xi
			v <- d %*% z$cov %*% t(d)                   # new VC matrix
			z$cov.t <- v
			z$se.t <- sqrt(diag(z$cov.t))               # new SEs
		}

		class(z) <- "gpd.fit"
		invisible(z)
	}

#------------------------------------------------------------------------------#
#              Algebraic observed information for GP(sigma, xi) fit             #
#------------------------------------------------------------------------------#

#' Observated information matrix for the Generalized Pareto distribution
#'
#' @param data data vector
#' @param scale scale parameter \eqn{sigma}
#' @param shape shape parameter \code{xi}
#' @param loc optional location parameter, corresponding to threshold
#'
#' @author Paul J. Northrop and Claire L. Coleman
#' @export
#' @return the observed information matrix
.gpd_obs_info <- function(data, scale, shape, loc=NULL) {
	y <- data
	if(!is.null(loc)){
		y <- y - loc
	}
	x <- shape; s <- scale
	i <- matrix(NA, 2, 2)
	i[1, 1] <- -sum((1 - (1 + x) * y * (2 * s + x * y) / (s + x * y)^2) /	s^2)
	i[1, 2] <- -sum(y * (1 - y / s) / (1 + x * y / s)^2 / s^2)
	i[2, 1] <- i[1, 2]
	i[2, 2] <-		sum(2 * log(1 + x * y / s) / x^3 - 2 * y / (s + x * y) / x^2 - (1 +
				1 / x) * y^2 / (s + x * y)^2)
	i
}
#Expected information based covariance matrix (Smith 1984)
gpd.vcov.mat <- function(data, scale, shape, loc=NULL){
		y <- data
	if(!is.null(loc)){
		y <- y - loc
	}
	info <- matrix(NA, 2, 2)
	info[1, 1] <- 1/scale^2
	info[1, 2] <- exp(-log(scale)-log(1+shape))
	info[2, 1] <- info[1, 2]
	info[2, 2] <- 2/(1+shape)
	info*scale^2*(1+shape)^2
}



#------------------------------------------------------------------------------#
#               Fit GP(sigma, xi) distribution using 2D optimisation           #
#------------------------------------------------------------------------------#

# Stolen from ismev; gradient function added.
#' @export
.gpd_2D_fit <-	function (xdat, threshold, npy = 365, ydat = NULL, sigl = NULL,
		shl = NULL, siglink = identity, shlink = identity, siginit = NULL,
		shinit = NULL, show = TRUE, method = "Nelder-Mead", maxit = 10000, do.fscale =
			FALSE, do.pscale = FALSE, ...) {
		z <- list()
		npsc <- length(sigl) + 1
		npsh <- length(shl) + 1
		n <- length(xdat)
		z$trans <- FALSE
		if (is.function(threshold))
			stop("`threshold' cannot be a function")
		u <- rep(threshold, length.out = n)
		if (length(unique(u)) > 1)
			z$trans <- TRUE
		xdatu <- xdat[xdat > u]
		xind <- (1:n)[xdat > u]
		u <- u[xind]
		in2 <- sqrt(6 * var(xdat)) / pi
		in1 <- mean(xdat, na.rm = TRUE) - 0.57722 * in2
		if (is.null(sigl)) {
			sigmat <- as.matrix(rep(1, length(xdatu)))
			if (is.null(siginit))
				siginit <- in2
		}
		else {
			z$trans <- TRUE
			sigmat <- cbind(rep(1, length(xdatu)), ydat[xind, sigl])
			if (is.null(siginit))
				siginit <- c(in2, rep(0, length(sigl)))
		}
		if (is.null(shl)) {
			shmat <- as.matrix(rep(1, length(xdatu)))
			if (is.null(shinit))
				shinit <- 0.1
		}
		else {
			z$trans <- TRUE
			shmat <- cbind(rep(1, length(xdatu)), ydat[xind, shl])
			if (is.null(shinit))
				shinit <- c(0.1, rep(0, length(shl)))
		}
		init <- c(siginit, shinit)
		z$model <- list(sigl, shl)
		z$link <- deparse(substitute(c(siglink, shlink)))
		z$threshold <- threshold
		z$nexc <- length(xdatu)
		z$data <- xdatu
		gpd.lik <- function(a) {
			sc <- siglink(sigmat %*% (a[seq(1, length = npsc)]))
			xi <- shlink(shmat %*% (a[seq(npsc + 1, length = npsh)]))
			y <- (xdatu - u) / sc
			y <- 1 + xi * y
			if (min(sc) <= 0){
				l <- 10^6
			} else {
				if (min(y) <= 0){
					l <- 10^6
				} else {
					l <- sum(log(sc)) + sum(log(y) * (1 / xi + 1))
				}
			}
			l
		}
		#
		gp.grad <- function(a) {
			# gradient of negated log-likelihood
			sigma <- a[1]; xi <- a[2]
			y <- xdatu - u
			n <- length(y)
			yy <- 1 + xi * y / sigma
						s1 <- -n * sigma^(-1) + (1 + xi) * sigma^(-2) * sum(y / yy)
			if(any(yy<0)){ #TODO FIX WHETHER THIS MAKES SENSE
				-c(s1,1e30)
				} else{
			s2 <- xi^(-2) * sum(log(yy)) - (1 + 1 / xi) * sigma^(-1) * sum(y / yy)
			- c(s1, s2)
				}
		}
		#
		x <- try(optim(init, gpd.lik, gr = gp.grad, hessian = TRUE, method = method,
				control = list(maxit = maxit, ...)
			)
		)
		if(is.character(x)){
		  z$conv <- 50
		  return(z)
		}
		sc <- siglink(sigmat %*% (x$par[seq(1, length = npsc)]))
		xi <- shlink(shmat %*% (x$par[seq(npsc + 1, length = npsh)]))
		z$conv <- x$convergence
		z$counts <- x$counts
		z$nllh <- x$value
		z$vals <- cbind(sc, xi, u)
		if (z$trans) {
			z$data <- -log(as.vector((1 + (xi * (
				xdatu - u
			)) / sc)^(-1 / xi)))
		}
		z$mle <- x$par
		z$rate <- length(xdatu) / n
		z$cov <- tryCatch(solve(x$hessian), error = function(e){matrix(NA, ncol(x$hessian),ncol(x$hessian))}) #TODO fix this
		suppressWarnings(z$se <- sqrt(diag(z$cov)))
		z$n <- n
		z$npy <- npy
		z$xdata <- xdat
		if (show) {
			if (z$trans)
				print(z[c(2, 3)])
			if (length(z[[4]]) == 1)
				print(z[4])
			print(z[c(5, 7)])
			if (!z$conv)
				print(z[c(8, 10, 11, 13)])
		}
		class(z) <- "gpd.fit"
		invisible(z)
	}

#------------------------------------------------------------------------------#
#       Fit GP(sigma, xi) distribution using 1D optimisation, using nlm        #
# ... to force gradients to be very close to zero at MLE                       #
#------------------------------------------------------------------------------#

.GP_1D_fit_nlm <-	function(xdat, threshold, init.val = NULL, calc.se = F, ...) {
		z.fit <- list()
		ex.data <- xdat[xdat > threshold] - threshold
		#
		GP.1D.negloglik <- function(phi) {
			# negated 1D loglikelihood
			k <- length(ex.data)
			if (phi == 0) {
				sigma.mle <- mean(ex.data)
				return(k * log(sigma.mle) + (1 / sigma.mle) * sum(ex.data))
			}
			zz <- 1 + phi * ex.data
			if (min(zz) <= 0)
				return(1e30)
			xi.of.phi <- mean(log(zz))
			sigma.of.phi <- xi.of.phi / phi
			if (sigma.of.phi <= 0)
				return(1e30)
			neg.log.lik <- 	k * log(sigma.of.phi) + (1 + 1 / xi.of.phi) * sum(log(1 + phi * ex.data))
			#
			attr(neg.log.lik, "gradient") <-	-k / phi + (1 + 1 / xi.of.phi) * sum(ex.data / zz)
			#
			attr(neg.log.lik, "hessian") <- 	k / phi^2 - xi.of.phi^(-2) * (sum(ex.data / zz))^2 / k - (1 + 1 /
						xi.of.phi) * sum(ex.data^2 / zz^2)
			#
			neg.log.lik
		}#................................# end of GP.1D.negloglik()
		#
		init <- ifelse(is.null(init.val), (2^0.1 - 1) / median(ex.data), init.val)
		f.scale <- GP.1D.negloglik(init)
		typ.size <- init
		suppressWarnings(res <- nlm(GP.1D.negloglik, init, fscale = f.scale, typsize = typ.size, iterlim=1000, check.analyticals =	FALSE, ...	))
		phi.hat <- res$estimate
		xi.hat <- mean(log(1 + res$estimate * ex.data))
		sigma.hat <- xi.hat / phi.hat
		z.fit$mle <- c(sigma.hat, xi.hat)
		z.fit$nllh <- res$minimum
		z.fit$gradient <- res$gradient
		z.fit$code <- res$code
		z.fit$counts["function"] <- res$iterations
		z.fit$convergence <- ifelse(res$code %in% 1:2, 0,1)
		if (calc.se)
			z.fit$se <-	sqrt(diag(solve(.gpd_obs_info(scale=sigma.hat, shape=xi.hat, data=ex.data))))
		invisible(z.fit)
	}


#------------------------------------------------------------------------------#
#                 GP fitting function of Grimshaw (1993)                       #
#------------------------------------------------------------------------------#

#' GP fitting function of Grimshaw (1993)
#'
#' Function for estimating parameters \code{k} and \code{a} for a random sample from a GPD.
#'
#' @author Paul J. Northrop and Claire L. Coleman
#' @param x sample values
#' @return a list with the maximum likelihood estimates of components \code{a} and \code{k}
.gpd_grimshaw <- function(x) {
	n  <- length(x)
	xq <- sort(x)
	xbar <- mean(x)
	sumx2 <- sum(x^2) / n
	x1 <- xq[1]
	xn <- xq[n]
	#  Find the local maxima/minima of the likelihood.
	epsilon <- 10^(-6) / xbar #  Initialize epsilon as the accuracy criterion
	#  The local maxima/minima must be found numerically by finding the zero(s) of h().
	#  Algorithm for finding the zero(s) of h().
	#  Any roots that exist must be within the interval (lobnd, hibnd).
	lobnd <- 2 * (x1 - xbar) / x1^2
	if (lobnd >= 0) {
		lobnd <- -epsilon
	}
	hibnd <- (1 / xn) - epsilon
	if (hibnd <= 0) {
		hibnd <- epsilon
	}
	#  If h''(0) > 0, look for one negative and one positive zero of h().
	#  If h''(0) < 0, look for two negative and two positive zeros of h().
	secderiv <- sumx2 - 2 * xbar^2  #{ Evaluate h''(0). }
	if (secderiv > 0) {
		#  Look for one negative and one positive zero of h().
		thzeros <- cbind(c(0, 0), c(0, 0))
		nzeros <- 2
		#  Begin with the initial value at lobnd.
		hlo <- (1 + sum(log(1 - lobnd * x)) / n) * (sum(1 / (1 - lobnd * x)) / n) - 1
		if (hlo < 0) {
			thlo <- lobnd       #{  Orient the search so h(thlo)<0  }
			thhi <- -epsilon
		} else {
			thlo <- -epsilon
			thhi <- lobnd
		}
		thzero <-		lobnd    #{  Initial value for modified Newton-Raphson is lobnd. }
		dxold <- abs(thhi - thlo)
		dx <- dxold
		temp1 <- sum(log(1 - thzero * x)) / n
		temp2 <- sum(1 / (1 - thzero * x)) / n
		temp3 <- sum(1 / (1 - thzero * x)^2) / n
		h <- (1 + temp1) * (temp2) - 1
		hprime <- (temp3 - temp2^2 - temp1 * (temp2 - temp3)) / thzero

		#  Newton-Raphson Algorithm to find the zero of the function h() for a given initial starting point.
		j <- 1
		maxiter <-	100  #{Maximum number of mod. Newton-Raphson iterations}
		while (j <= maxiter) {
			#  Determine whether it is better to use Bisection (if N-R is
			#  out of range or not decreasing fast enough) or Newton-Raphson.
			c1 <- (((thzero - thhi) * hprime - h) * ((thzero - thlo) * hprime - h) >=	0)
			c2 <- (abs(2 * h) > abs(dxold * hprime))
			if (c1 + c2 >= 1) {
				dxold <- dx
				dx <- (thhi - thlo) / 2
				thzero <- thlo + dx
				if (thlo == thzero) {
					#{Change in root is negligible}
					j <- 1000
				}
			} else{
				dxold <- dx
				dx <- h / hprime
				temp <- thzero
				thzero <- thzero - dx
				if (temp == thzero) {
					#{Change in root is negligible}
					j <- 1001
				}
			}
			#  Determine if convergence criterion is met.
			if (abs(dx) < epsilon * abs(thlo + thhi) / 2) {
				j <- 999
			}
			temp1 <- sum(log(1 - thzero * x)) / n
			temp2 <- sum(1 / (1 - thzero * x)) / n
			temp3 <- sum(1 / (1 - thzero * x)^2) / n
			h <- (1 + temp1) * (temp2) - 1
			hprime <- (temp3 - temp2^2 - temp1 * (temp2 - temp3)) / thzero
			if (h < 0) {
				#{Maintain the bracket on the root}
				thlo <- thzero
			}  else{
				thhi <- thzero
			}
			j <- j + 1
		}

		if (j > maxiter + 1) {
			thzeros[1, ] <- cbind(thzero, j)
		}
		#  Begin with the initial value at hibnd.
		hlo <- (1 + sum(log(1 - epsilon * x)) / n) * (sum(1 / (1 - epsilon * x)) /	n) - 1
		if (hlo < 0) {
			thlo <- epsilon       #{  Orient the search so h(thlo)<0  }
			thhi <- hibnd
		} else{
			thlo <- hibnd
			thhi <- epsilon
		}
		thzero <-	hibnd    #{  Initial value for modified Newton-Raphson is hibnd. }
		dxold <- abs(thhi - thlo)
		dx <- dxold
		temp1 <- sum(log(1 - thzero * x)) / n
		temp2 <- sum(1 / (1 - thzero * x)) / n
		temp3 <- sum(1 / (1 - thzero * x)^2) / n
		h <- (1 + temp1) * (temp2) - 1
		hprime <- (temp3 - temp2^2 - temp1 * (temp2 - temp3)) / thzero

		#  Newton-Raphson Algorithm to find the zero of the function h()
		#  for a given initial starting point.
		j <- 1
		maxiter <- 100  #{Maximum number of mod. Newton-Raphson iterations}
		while (j <= maxiter) {
			#  Determine whether it is better to use Bisection (if N-R is
			#  out of range or not decreasing fast enough) or Newton-Raphson.
			c1 <- (((thzero - thhi) * hprime - h) * ((thzero - thlo) * hprime - h) >=	0)
			c2 <- (abs(2 * h) > abs(dxold * hprime))
			if (c1 + c2 >= 1) {
				dxold <- dx
				dx <- (thhi - thlo) / 2
				thzero <- thlo + dx
				if (thlo == thzero) {
					#{Change in root is negligible}
					j <- 1000
				}
			}  else{
				dxold <- dx
				dx <- h / hprime
				temp <- thzero
				thzero <- thzero - dx
				if (temp == thzero) {
					#{Change in root is negligible}
					j <- 1001
				}
			}
			#  Determine if convergence criterion is met.
			if (abs(dx) < epsilon * abs(thlo + thhi) / 2) {
				j <- 999
			}
			temp1 <- sum(log(1 - thzero * x)) / n
			temp2 <- sum(1 / (1 - thzero * x)) / n
			temp3 <- sum(1 / (1 - thzero * x)^2) / n
			h <- (1 + temp1) * (temp2) - 1
			hprime <- (temp3 - temp2^2 - temp1 * (temp2 - temp3)) / thzero
			if (h < 0) {
				#{Maintain the bracket on the root}
				thlo <- thzero
			} else{
				thhi <- thzero
			}
			j <- j + 1
		}
		if (j > maxiter + 1) {
			thzeros[2, ] <- cbind(thzero, j)
		}
	}  else{ #if
		##  Look for two negative and two positive zeros of h().
		thzeros <- matrix(rep(0, 8), ncol = 2)
		nzeros <- 4
		#  Begin with the initial value at lobnd.
		hlo <- (1 + sum(log(1 - lobnd * x)) / n) * (sum(1 / (1 - lobnd * x)) /
				n) - 1
		if (hlo < 0) {
			thlo <- lobnd       #{  Orient the search so h(thlo)<0  }
			thhi <- -epsilon
		}  else{
			thlo <- -epsilon
			thhi <- lobnd
		}
		thzero <-	lobnd    #{  Initial value for modified Newton-Raphson is lobnd. }
		dxold <- abs(thhi - thlo)
		dx <- dxold
		temp1 <- sum(log(1 - thzero * x)) / n
		temp2 <- sum(1 / (1 - thzero * x)) / n
		temp3 <- sum(1 / (1 - thzero * x)^2) / n
		h <- (1 + temp1) * (temp2) - 1
		hprime <- (temp3 - temp2^2 - temp1 * (temp2 - temp3)) / thzero

		##  Newton-Raphson Algorithm to find the zero of the function h() for a given initial starting point.
		j <- 1
		maxiter <- 100  #{Maximum number of mod. Newton-Raphson iterations}
		while (j <= maxiter) {
			#  Determine whether it is better to use Bisection (if N-R is
			#  out of range or not decreasing fast enough) or Newton-Raphson.
			c1 <- (((thzero - thhi) * hprime - h) * ((thzero - thlo) * hprime - h) >=	0)
			c2 <- (abs(2 * h) > abs(dxold * hprime))
			if (c1 + c2 >= 1) {
				dxold <- dx
				dx <- (thhi - thlo) / 2
				thzero <- thlo + dx
				if (thlo == thzero) {
					#{Change in root is negligible}
					j <- 1000
				}
			}  else{
				dxold <- dx
				dx <- h / hprime
				temp <- thzero
				thzero <- thzero - dx
				if (temp == thzero) {
					#{Change in root is negligible}
					j <- 1001
				}
			}
			#
			#  Determine if convergence criterion is met.
			#
			if (abs(dx) < epsilon * abs(thlo + thhi) / 2) {
				j <- 999
			}
			temp1 <- sum(log(1 - thzero * x)) / n
			temp2 <- sum(1 / (1 - thzero * x)) / n
			temp3 <- sum(1 / (1 - thzero * x)^2) / n
			h <- (1 + temp1) * (temp2) - 1
			hprime <- (temp3 - temp2^2 - temp1 * (temp2 - temp3)) / thzero
			if (h < 0) {
				#{Maintain the bracket on the root}
				thlo <- thzero
			} else{
				thhi <- thzero
			}
			j <- j + 1
		}
		if (j > maxiter + 1) {
			thzeros[1, ] <- cbind(thzero, j)
		}
		#  Look at the derivative to determine where the second root lies.
		#   If h'(0)>0, second root lies between thzero and -epsilon.
		#   If h'(0)<0, second root lies between lobnd and thzero.
		temp1 <- sum(log(1 - thzero * x)) / n
		temp2 <- sum(1 / (1 - thzero * x)) / n
		temp3 <- sum(1 / (1 - thzero * x)^2) / n
		hprime <- (temp3 - temp2^2 - temp1 * (temp2 - temp3)) / thzero
		if (hprime > 0) {
			#  h'(0)>0, so the second zero lies between thzero and -epsilon.
			#  Establish Initial Values.
			thlo <- thzero
			thhi <- -epsilon
			thzero <- thhi
			dx <- thlo - thhi
			j <- 1
			maxiter <- 100  #{Maximum number of bisection iterations}
			while (j <= maxiter) {
				dx <- .5 * dx
				thmid <- thzero + dx
				hmid <- (1 + sum(log(1 - thmid * x)) / n) * (sum(1 / (1 - thmid * x)) /	n) - 1
				if (hmid < 0) {
					thzero <- thmid
				} else if (hmid == 0) {
					#{Zero of h() has been found}
					j <- 999
				}
				if (abs(dx) < epsilon * abs(thlo + thhi) / 2) {
					j <- 999
				}
				j <- j + 1
			}

			if (j > maxiter + 1) {
				thzeros[2, ] <- cbind(thzero, j)
			}
		}  else{
			#  h'(0)<0, so the second zero lies between lobnd and thzero.
			#  Establish Initial Values.
			#
			thlo  <- lobnd
			thhi <- thzero
			thzero <- thlo
			dx <- thhi - thlo
			j <- 1
			maxiter <- 100  #{Maximum number of bisection iterations}
			while (j <= maxiter) {
				dx <- .5 * dx
				thmid <- thzero + dx
				hmid <- (1 + sum(log(1 - thmid * x)) / n) * (sum(1 / (1 - thmid * x)) / n) - 1
				if (hmid < 0) {
					thzero <- thmid
				} else if (hmid == 0) {
					#{Zero of h() has been found}
					j <- 999
				}
				if (abs(dx) < epsilon * abs(thlo + thhi) / 2) {
					j <- 999
				}
				j <- j + 1
			}

			if (j > maxiter + 1) {
				thzeros[2, ] <- cbind(thzero, j)
			}
		}
		#  Begin with the initial value at hibnd.
		hlo <- (1 + sum(log(1 - epsilon * x)) / n) * (sum(1 / (1 - epsilon * x)) /
				n) - 1
		if (hlo < 0) {
			thlo <- epsilon       #{  Orient the search so h(thlo)<0  }
			thhi <- hibnd
		} else{
			thlo <- hibnd
			thhi <- epsilon
		}
		thzero <-
			hibnd    #{  Initial value for modified Newton-Raphson is hibnd. }
		dxold <- abs(thhi - thlo)
		dx <- dxold
		temp1 <- sum(log(1 - thzero * x)) / n
		temp2 <- sum(1 / (1 - thzero * x)) / n
		temp3 <- sum(1 / (1 - thzero * x)^2) / n
		h <- (1 + temp1) * (temp2) - 1
		hprime <- (temp3 - temp2^2 - temp1 * (temp2 - temp3)) / thzero
		#  Newton-Raphson Algorithm to find the zero of the function h()
		#  for a given initial starting point.
		j <- 1
		maxiter <- 100  #{Maximum number of mod. Newton-Raphson iterations}
		while (j <= maxiter) {
			#  Determine whether it is better to use Bisection (if N-R is
			#  out of range or not decreasing fast enough) or Newton-Raphson.
			c1 <- (((thzero - thhi) * hprime - h) * ((thzero - thlo) * hprime - h) >=	0)
			c2 <- (abs(2 * h) > abs(dxold * hprime))
			if (c1 + c2 >= 1) {
				dxold <- dx
				dx <- (thhi - thlo) / 2
				thzero <- thlo + dx
				if (thlo == thzero) {
					#{Change in root is negligible}
					j <- 1000
				}
			} else{
				dxold <- dx
				dx <- h / hprime
				temp <- thzero
				thzero <- thzero - dx
				if (temp == thzero) {
					#{Change in root is negligible}
					j <- 1001
				}
			}
			#  Determine if convergence criterion is met.
			if (abs(dx) < epsilon * abs(thlo + thhi) / 2) {
				j <- 999
			}
			temp1 <- sum(log(1 - thzero * x)) / n
			temp2 <- sum(1 / (1 - thzero * x)) / n
			temp3 <- sum(1 / (1 - thzero * x)^2) / n
			h <- (1 + temp1) * (temp2) - 1
			hprime <- (temp3 - temp2^2 - temp1 * (temp2 - temp3)) / thzero
			if (h < 0) {
				#{Maintain the bracket on the root}
				thlo <- thzero
			} else{
				thhi <- thzero
			}
			j <- j + 1
		}

		if (j > maxiter + 1) {
			thzeros[3, ] <- cbind(thzero, j)
		}
		#  Look at the derivative to determine where the second root lies.
		#   If h'(0)>0, second root lies between thzero and hibnd.
		#   If h'(0)<0, second root lies between epsilon and thzero.
		temp1 <- sum(log(1 - thzero * x)) / n
		temp2 <- sum(1 / (1 - thzero * x)) / n
		temp3 <- sum(1 / (1 - thzero * x)^2) / n
		hprime <- (temp3 - temp2^2 - temp1 * (temp2 - temp3)) / thzero
		if (hprime > 0) {
			#  h'(0)>0, so the second zero lies between thzero and hibnd.
			#  Establish Initial Values.
			thlo <- thzero
			thhi <- hibnd
			thzero <- thhi
			dx <- thlo - thhi

			j <- 1
			maxiter <- 100  #{Maximum number of bisection iterations}
			while (j <= maxiter) {
				dx <- .5 * dx
				thmid <- thzero + dx
				hmid <- (1 + sum(log(1 - thmid * x)) / n) * (sum(1 / (1 - thmid * x)) /	n) - 1
				if (hmid < 0) {
					thzero <- thmid
				} else if (hmid == 0) {
					#{Zero of h() has been found}
					j <- 999
				}
				if (abs(dx) < epsilon * abs(thlo + thhi) / 2) {
					j <- 999
				}
				j <- j + 1
			}

			if (j > maxiter + 1) {
				thzeros[4, ] <- cbind(thzero, j)
			}
		} else{
			#  h'(0)<0, so the second zero lies between epsilon and thzero.
			#  Establish Initial Values.
			thlo <- epsilon
			thhi <- thzero
			thzero <- thlo
			dx <- thhi - thlo

			j <- 1
			maxiter <- 100  #{Maximum number of bisection iterations}
			while (j <= maxiter) {
				dx <- .5 * dx
				thmid <- thzero + dx
				hmid <- (1 + sum(log(1 - thmid * x)) / n) * (sum(1 / (1 - thmid * x)) / n) - 1
				if (hmid < 0) {
					thzero <- thmid
				} else if (hmid == 0) {
					#{Zero of h() has been found}
					j <- 999
				}
				if (abs(dx) < epsilon * abs(thlo + thhi) / 2) {
					j <- 999
				}
				j <- j + 1
			}

			if (j > maxiter + 1) {
				thzeros[4, ] <- cbind(thzero, j)
			}
		}
	}
	#  Of the candidate zero(s) of h(), determine whether they correspond
	#  to a local maximum or minimum of the log-likelihood.
	#  Eliminate any non-convergent roots}
	thetas <- thzeros[thzeros[, 2] > maxiter + 1, ]
	nzeros <- nrow(thetas)
	proll <- rep(0, nzeros)
	mles <- matrix(rep(0, 4 * nzeros), ncol = 4)
	i <- 1
	while (i <= nzeros) {
		temp1 <- sum(log(1 - thetas[i, 1] * x))
		mles[i, 1] <- -temp1 / n
		mles[i, 2] <- mles[i, 1] / thetas[i, 1]
		mles[i, 3] <- -n * log(mles[i, 2]) + (1 / mles[i, 1] - 1) * temp1
		mles[i, 4] <- 999
		i <- i + 1
	}
	ind <- 1:length(mles[, 4])
	ind <- ind[mles[, 4] == 999]
	if (sum(ind) == 0) {
		#{ Check to see if there are any local maxima. }
		nomle <- 0          #{ If not, there is no mle. }
	} else{
		nomle <- 1
	}
	if (nomle != 0) {
		mles  <- mles[ind, ]
		nmles <- nrow(mles)
		#  Add the boundary value where k=1 to the candidates for the
		#  maximum of the log-likelihood.
		mles <- rbind(mles, c(1, xn, -n * log(xn), 999))
		nmles <- nmles + 1
		#  Choose of the candidate mles whichever has the largest log-likelihood.
		maxlogl <- max(mles[, 3])
		ind <- order(mles[, 3])
		ind <- ind[nmles]
		k <- mles[ind, 1]
		#  label(k, 'GPD mle of k')
		a <- mles[ind, 2]
		conv <- 0
		#  label(a, 'GPD mle of a')
	} else{
		conv <- 50
		#  No Maximum Likelihood Estimators were found.
		k <- NA
		a <- NA
	}
	list(k = k, a = a, conv=conv )
}

".Zhang_Stephens_posterior" <-
function(x) {# x: sample data from the GPD
  n <- length(x);
  x <- sort(x)
  lx <- function(b, x){
    k <- -mean(log(1-b*x));
    log(b/k)+k-1
    }
  m <- 20 + floor(sqrt(n))
  b <- 1/x[n]+(1-sqrt(m/(1:m-.5)))/3/x[floor(n/4+.5)]
  L <- sapply(1:m, function(i){ n*lx(b[i],x)})
  w <- sapply(1:m, function(i){1/sum(exp(L-L[i]))})
  b <- sum(b*w);
  xi <- mean(log(1-b*x)); sigma <- -xi/b
  mle <- c(sigma, xi)
  names(mle) <- c("scale", "shape")
  list(mle = mle)
}

".Zhang_posterior" <- function(x) {
  n <- length(x) ;
  x <- sort(x)
  lx <- function(b, x) {
    k <- -mean(log(1-b*x))
    if (b==0){
      k-1-log(mean(x))
    } else{
      k-1+log(b/k)
    }
  }
  p <- (3:9)/10 ; xp <- x[round(n*(1-p)+.5)]
  m <- 20+round(n^.5) ;
  xq <- x[round(n*(1-p*p)+.5)]
  k <- log(xq/xp-1,p) ; a <- k*xp/(1-p^k)
  a[k==0] <- (-xp/log(p))[k==0] ;
  k <- -1
  b <- w <- L <- (n-1)/(n+1)/x[n]-(1-((1:m-.5)/m)^k)/k/median(a)/2
  L <- sapply(1:m, function(i) n*lx(b[i],x))
  w <- sapply(1:m, function(i) 1/sum(exp(L-L[i])))
  b <- sum(b*w)
  k <- -mean(log(1-b*x))
  mle <- c(k/b, -k)
  names(mle) <- c("scale", "shape")
  list(mle = mle)
}



#' Peaks-Over-Threshold Modelling using the Generalized Pareto Distribution
#'
#' Numerical optimization of the Generalized Pareto distribution over a
#' high threshold.
#'
#' @param xdat a numeric vector of data to be fitted.
#' @param threshold the chosen threshold.
#' @param show logical; if \code{TRUE} (the default), print details of the fit.
#' @param method the method to be used. See \bold{Details}. Can be abbreviated.
#' @param MCMC \code{NULL} for frequentist estimates, otherwise a boolean or a list with parameters passed. If \code{TRUE}, runs a Metropolis-Hastings sampler to get posterior mean estimates. Can be used to pass arguments \code{niter}, \code{burnin} and \code{thin} to the sampler as a list.
#' @seealso \code{\link[evd]{fpot}} and \code{\link[ismev]{gpd.fit}}
#'
#' @details The default method is \code{"Grimshaw"}, consisting in maximization of the profile likelihood for the scale.
#' Other options for maximization of the profile likelihood are \code{nlm} and \code{optim}, which use respectively \code{\link[stats]{nlm}} and \code{\link[stats]{optim}}. Method \code{"ismev"} is the two-dimensional optimization routine \code{\link[ismev]{gpd.fit}} from the \code{\link[ismev]{ismev}} library, with in addition the algebraic gradient.
#' The approximate Bayesian methods (\code{"zs"} and \code{"zhang"}) are extracted respectively from Zhang and Stephens (2009) and Zhang (2010) and consists of a approximate posterior mean calculated via importance
#' sampling assuming a GPD prior is placed on the parameter of the profile likelihood.
#' @note Some of the internal functions (which are hidden from the user) allow for modelling of the parameters using covariates. This is not currently implemented within \code{gp.fit}, but users can call internal functions should they wish to use these features.
#' @author Paul J. Northrop and Claire L. Coleman for the frequentist functions.
#' Zhang and Stephens (2009) and Zhang (2010) for the \code{zs} and \code{zhang} approximate methods and L. Belzile for the wrapper and MCMC samplers.
#'
#' @references Davison, A.C. (1984). Modelling excesses over high thresholds, with an application, in
#' \emph{Statistical extremes and applications}, J. Tiago de Oliveira (editor), D. Reidel Publishing Co., 461--482.
#' @references Grimshaw, S.D. (1993). Computing Maximum Likelihood Estimates for the Generalized
#'  Pareto Distribution, \emph{Technometrics}, \bold{35}(2), 185--191.
#' @references Northrop, P.J. and C. L. Coleman (2014). Improved threshold diagnostic plots for extreme value
#' analyses, \emph{Extremes}, \bold{17}(2), 289--303.
#' @references Zhang, J. (2010). Improving on estimation for the generalized Pareto distribution, \emph{Technometrics} \bold{52}(3), 335--339.
#' @references Zhang, J.  and M.A. Stephens (2009). A new and efficient estimation method for the generalized Pareto distribution.
#' \emph{Technometrics} \bold{51}(3), 316--325.
#'
#'
#' @return If \code{method} is neither \code{"zs"} nor \code{"zhang"}, a list containing the following components:
#' \itemize{
#' \item \code{estimate} a vector containing all parameters (optimized and fixed).
#' \item \code{std.err} a vector containing the standard errors.
#' \item \code{var.cov} the variance covariance matrix, obtained as the numerical inverse of the observed information matrix.
#' \item \code{threshold} the threshold.
#' \item \code{method} the method used to fit the parameter. See details.
#' \item \code{deviance} the deviance at the maximum likelihood estimates.
#' \item \code{nat} number of points lying above the threshold.
#' \item \code{pat} proportion of points lying above the threshold.
#' \item \code{convergence} components taken from the list returned by \code{\link[stats]{optim}}.
#' Values other than \code{0} indicate that the algorithm likely did not converge (in particular 1 and 50).
#' \item \code{counts} components taken from the list returned by \code{\link[stats]{optim}}.
#' }
#' Otherwise, a list containing
#' \itemize{
#' \item \code{threshold} the threshold.
#' \item \code{method} the method used to fit the parameter. See \bold{Details}.
#' \item \code{nat} number of points lying above the threshold.
#' \item \code{pat} proportion of points lying above the threshold.
#' \item \code{approx.mean} a vector containing containing the approximate posterior mean estimates.
#' }
#' and in addition if MCMC is neither \code{FALSE}, nor \code{NULL}
#' \itemize{
#' \item \code{post.mean} a vector containing the posterior mean estimates.
#' \item \code{post.se} a vector containing the posterior standard error estimates.
#' \item \code{accept.rate} proportion of points lying above the threshold.
#' \item \code{niter} length of resulting Markov Chain
#' \item \code{burnin} amount of discarded iterations at start, capped at 10000.
#' \item \code{thin} thinning integer parameter describing
#' }
#'
#' @export
#'
#' @examples
#' library(ismev)
#' data(rain)
#' threshold <- quantile(rain,0.9)
#' gp.fit(rain, threshold, method="Grimshaw")
#' gp.fit(rain, threshold, method="zs")
gp.fit <- function(xdat, threshold, method=c("Grimshaw","nlm","optim","ismev","zs","zhang"), show=TRUE, MCMC=NULL){
	xi.tol = 1e-3
	#Optimization of model, depending on routine
	method <- match.arg(method)
	if(!is.null(MCMC) && ! method %in% c("zs","zhang")) warning("Ignoring argument `MCMC` for frequentist estimation")
	if(missing(method)){
		method="Grimshaw"
	}
				if (method == "nlm") {
				temp <- .Zhang_Stephens_posterior(xdat)
# 				temp <-	.gpd_1D_fit(xdat, threshold, show = F, xi.tol = xi.tol, calc.se = FALSE)
				# threshold 1D max, algebraic Hessian
# 				#If algorithm failed to converge with initial values, switch to Grimshaw
# 				if(temp$conv != 0 || any(is.nan(temp$mle))){ #algorithm
# 					warning("Algorithm did not converge. Switching method to Grimshaw")
# 					#if( temp$mle[1] < 0 || temp$mle[2]>10 || temp$mle[2]<3 TODO fix
# 				  temp <- .gpd_grimshaw(xdat[xdat > threshold] - threshold)
# 				  temp$mle <- c(temp$a, -temp$k)
# 				  #If values are numeric, but the algorithm diverged
# 				} else if (is.numeric(temp$mle) && !is.nan(temp$mle)){
# 					if( temp$mle[1] < 0 || temp$mle[2]>10 || temp$mle[2]<3){
# 						temp <- .gpd_grimshaw(xdat[xdat > threshold] - threshold)
# 				  temp$mle <- c(temp$a, -temp$k)
# 					}
# 				}
				temp <-	.GP_1D_fit_nlm(xdat, threshold, init.val = temp$mle[2] / temp$mle[1],
				                       gradtol = 1e-10, steptol = 1e-5, calc.se = FALSE)
        ifelse(temp$code <2, temp$conv <- 0, temp$conv <- 50)
        # 1D max, use nlm to get gradients v close to zero
			} else if (method == "ismev") {
				temp <-	.gpd_2D_fit(xdat, threshold, show = FALSE) # ismev GP fitting function
				if(temp$conv != 0){#algorithm failed to converge
				  warning("Algorithm did not converge. Switching method to Grimshaw")
				  temp <- .gpd_grimshaw(xdat[xdat > threshold] - threshold)
				  temp$mle <- c(temp$a, -temp$k)
				}
				temp <-	.gpd_2D_fit(xdat, threshold, show = FALSE, siginit = temp$mle[1], shinit = temp$mle[2],
				                    method = "BFGS", reltol = 1e-30, abstol = 1e-30)
			} else if (method == "optim") {
				temp <-	.gpd_1D_fit(xdat, threshold, show = FALSE, xi.tol = xi.tol) # 1D max, algebraic Hessian
				if(temp$conv != 0){#algorithm failed to converge
				  warning("Algorithm did not converge. Switching method to Grimshaw")
				  temp <- .gpd_grimshaw(xdat[xdat > threshold] - threshold)
				  temp$mle <- c(temp$a, -temp$k)
				}
				temp <-	.gpd_1D_fit(xdat, threshold, show = FALSE, xi.tol = xi.tol, phi.input = temp$mle[2] / temp$mle[1],
				                    reltol = 1e-30, abstol = 1e-30) #1D max, algebraic Hessian
			} else if (method=="zs" || method=="zhang") {
			  xdat = xdat[xdat > threshold]-threshold
				temp <- switch(method,
				            zs = .Zhang_Stephens_posterior(xdat),
				            zhang = .Zhang_posterior(xdat)
				)
				if(!is.null(MCMC) && MCMC!=FALSE){
				if(is.logical(MCMC) && !is.na(MCMC) && MCMC){ #if TRUE
						burn = 2000; thin = 1; niter= 10000
					} else if(is.list(MCMC)){
						burn <- ifelse(is.null(MCMC$burnin),2000,MCMC$burnin)
						thin <- ifelse(is.null(MCMC$thin),1,MCMC$thin)
						niter <- ifelse(is.null(MCMC$niter),10000,MCMC$niter)
					}
				bayespost <- switch(method,
				                  zs = Zhang_Stephens(xdat, init=-temp$mle[2]/temp$mle[1], burnin=burn, thin=thin,
				                  										niter=niter, method=1),
				                  zhang = Zhang_Stephens(xdat, init=-temp$mle[2]/temp$mle[1], burnin=burn,
				                  											 thin=thin, niter=niter, method=2)
				)
				#Copying output for formatting and printing
				post.mean <- bayespost$summary[1,]
				post.se <- sqrt(bayespost$summary[2,])
				names(post.mean) <- names(post.se) <- c("scale","shape")
				post <- structure(list(
					method = method,
					threshold = threshold,
					nat = sum(xdat > threshold),
					pat=sum(xdat > threshold) / length(xdat),
					approx.mean = temp$mle,
					post.mean = post.mean,
					post.se = post.se,
					accept.rate = bayespost$rate,
					thin = bayespost$thin,
					burnin=bayespost$burnin,
					niter=bayespost$niter), class="gpdbayes")

				} else{
					post <- structure(list(
						method = method,
						threshold = threshold,
						nat = sum(xdat > threshold),
						pat=sum(xdat > threshold) / length(xdat),
						approx.mean = temp$mle), class="gpdbayes")

				}
				if(show) print(post)
			return(invisible(post))
     }
		if (method == "Grimshaw") {
				yy <- xdat[xdat > threshold] - threshold # thresholds excesses
				pjn <-	.gpd_grimshaw(yy)  # Grimshaw (1993) function, note: k is -xi, a is sigma
				temp <- list()
				temp$mle <- c(pjn$a, -pjn$k)        # mle for (sigma, xi)
				sc <- rep(temp$mle[1], length(yy));
				xi <- temp$mle[2]
				temp$nllh <- sum(log(sc)) + sum(log(1 + xi * yy / sc) * (1 / xi + 1))
				temp$conv <- pjn$conv
			}


		#Collecting observations from temp and formatting the output
	invobsinfomat <- tryCatch(solve(.gpd_obs_info(data=xdat[xdat>threshold], scale=temp$mle[1],
																shape=temp$mle[2], loc=threshold)), error= function(e){"notinvert"}, warning= function(w) w)

		if(invobsinfomat =="notinvert" || all(is.nan(invobsinfomat))){
			warning("Cannot calculate standard error based on observed information")
			if(!is.null(temp$se)){std.errors <- diag(temp$se)} else{std.errors <- diag(rep(NA,2))}
		} else if(!is.null(temp$mle) && temp$mle[2] > -0.5 && temp$conv==0){#If the MLE was returned
		std.errors <- sqrt(diag(invobsinfomat))
		} else{
			warning("Cannot calculate standard error based on observed information")
			std.errors <- rep(NaN, 2)
		}
		if(temp$mle[2] < -1 && temp$conv == 0){
			warning("The MLE is not a solution to the score equation for `xi < -1'")
		}
		names(temp$mle) <- names(std.errors) <- c("scale","shape")
		output <- structure(list(threshold = threshold,
				estimate = temp$mle,
				std.err = std.errors,
				var.cov = invobsinfomat,
				threshold = threshold,
				method = method,
				deviance = 2*temp$nllh,
				nat = sum(xdat > threshold),
				pat=sum(xdat > threshold) / length(xdat),
				convergence = temp$conv,
				counts = temp$counts
			), class = "gpd")
		if(show){
			print(output)
		}
			invisible(output)
}

# #' @param x A fitted object of class \code{gpd}.
# #' @param digits Number of digits to display in \code{print} call.
# #' @param ... Additional argument passed to \code{print}.
# #' @rdname gp.fit
# #' @S3method print gpd
# #' @exportMethods
#' @export
"print.gpd" <-  function(x, digits = max(3, getOption("digits") - 3), ...){
	 	cat("Method:", x$method, "\n")
    cat("Deviance:", x$deviance, "\n")

    cat("\nThreshold:", round(x$threshold, digits), "\n")
    cat("Number Above:", x$nat, "\n")
    cat("Proportion Above:", round(x$pat, digits), "\n")

    cat("\nEstimates\n")
    print.default(format(x$estimate, digits = digits), print.gap = 2,
        quote = FALSE,...)
    if(!is.null(x$std.err) && x$estimate[1] > -0.5) {
    cat("\nStandard Errors\n")
    print.default(format(x$std.err, digits = digits), print.gap = 2,
        quote = FALSE,...)
    }
    cat("\nOptimization Information\n")
    cat("  Convergence:", x$convergence, "\n")
    if(x$method!="Grimshaw"){
    cat("  Function Evaluations:", x$counts["function"], "\n")
    if(!is.null(x$counts["gradient"]) && x$method!="nlm")
        cat("  Gradient Evaluations:", x$counts["gradient"], "\n")
    cat("\n")
    }
    invisible(x)
}

# #' @param x A fitted object of class \code{gpdbayes}.
# #' @param digits Number of digits to display in \code{print} call.
# #' @param ... Additional argument passed to \code{print}.
# #' @rdname gp.fit
# #' @S3method print gpdbayes
# #' @exportMethods
#' @export
"print.gpdbayes" <-  function(x, digits = max(3, getOption("digits") - 3), ...){
	cat("\nMethod:", switch(x$method,zs="Zhang and Stephens",zhang="Zhang"), "\n")
	cat("\nThreshold:", round(x$threshold, digits), "\n")
	cat("Number Above:", x$nat, "\n")
	cat("Proportion Above:", round(x$pat, digits), "\n")

cat("\nApproximate posterior mean estimates\n")
print.default(format(x$approx.mean, digits = 3), print.gap = 2, quote = FALSE)
if(!is.null(x$post.mean)){
	cat("\nPosterior mean estimates\n")
	print.default(format(x$post.mean, digits = 3), print.gap = 2, quote = FALSE)
	cat("\nMonte Carlo standard errors\n")
	print.default(format(x$post.se, digits = 3), print.gap = 2, quote = FALSE)
	cat("\nEstimates based on an adaptive MCMC\n Runs:   ",x$niter,
			"\n Burnin: ", x$burnin,"\n Acceptance rate:",
			round(x$accept.rate,digits=2),"\n Thinning:",x$thin, "\n")

}
}
