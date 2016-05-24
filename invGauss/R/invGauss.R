invGauss <- function(formula.mu, formula.c = ~1, link.mu = identity, 
 link.c = exp, data, mu = TRUE, tau = TRUE, verbose = FALSE, protect = FALSE, 
 opti.method = "BFGS", use.gradient = TRUE, link.mu.deriv = function(x) 1, 
 link.c.deriv = exp)
{
#
## PREPARE
#if(!is.element(opti.method, c("nlmin", "nlminb", "ms", "optim"))) stop("invalid opti.method!")
if(length(opti.method) > 1) stop("Cannot run several methods at once")
.mmu <- missing(link.mu)
.mmud <- missing(link.mu.deriv)
.mc <- missing(link.c)
.mcd <- missing(link.c.deriv)
if(use.gradient & (
		((!.mmu & .mmud) | (.mmu & !.mmud)) | 
		((!.mc & .mcd) | (.mc & !.mcd))
   )) stop("both the link function and its derivative should be supplied when use.gradient = TRUE")
#
## EXTRACT time, status AND covariates
model.mu <- model.frame(formula.mu, data = data)
response <- model.extract(model.mu, "response")
.time <- response[, "time"]
status <- response[, "status"]
#
## test if censoring
no.censor <- all(status == 1)
#
## test if events
if(all(status != 1)) stop("No events!")	#
#
## prepare covariate matrices
covar.c <- model.matrix(formula.c, data = data)
covar.mu <- model.matrix(formula.mu, data = data)
#
ncovar.c <- ncol(covar.c)
ncovar.mu <- ncol(covar.mu)
#
## names for coefficients
termnames <- c("tau", paste(dimnames(covar.mu)[[2]], ".mu", sep = ""), paste(dimnames(covar.c)[[2]], ".c", sep = ""))	# xxx
termnames <- gsub("(Intercept)", "Intercept", termnames, fixed = T)
termnames <- make.names(termnames)  # should (usually) make them useable for data.frames
#
## MAKING AVAILABLE VARIABLES SPLIT BY status
covar.c0 <- covar.c[status == 0,  , drop = F]
covar.c1 <- covar.c[status == 1,  , drop = F]
.time0 <- .time[status == 0]
.time1 <- .time[status == 1]
covar.c.10 <- rbind(covar.c1, covar.c0)
#
covar.mu0 <- covar.mu[status == 0,  , drop = F]
covar.mu1 <- covar.mu[status == 1,  , drop = F]
covar.mu.10 <- rbind(covar.mu1, covar.mu0)

use.hess <- F
if(use.hess){
#			# ikke klargjort for hess over mu.. burde ikke vaere saa vanskelig
#			#
#			temp <- covar.c.10[, c(1, 1, 1:ncovar.c)]
#			covar.c.10.hessian <- array(NA, dim = dim(temp)[c(1, 2, 2)])
#			for(j in 1:dim(covar.c)[1]) covar.c.10.hessian[j,  ,  ] <- temp[j,  ] %*% t(temp[j,  #			  ])
}
#
#
i <- 0
.f.bB <- function(x, use.grad){
	##
	## COMPUTE CONTRIBUTIONS FROM f.fb and f.B TO LOG-LIKELIHOOD
	#
	## compute for events
	if(verbose) print(x, digits = 2)	#
	mu1.link <- link.mu(covar.mu1 %*% x[2:(1 + ncovar.mu)]) #
	c1.link <- link.c(covar.c1 %*% x[(2 + ncovar.mu):(1 + ncovar.mu + ncovar.c)])
	if(protect)
		c1.link <- pmax(c1.link, 0)	#
	if(!use.grad){
		o.fb <- f.fb(t = .time1, mu = mu1.link, c1 = c1.link, tau = x[1])
	}else{
		o.fb <- f.fb.neglog.deriv(t = .time1, mu = mu1.link, c1 = c1.link, tau = x[1])	
	}
	.test <- any(o.fb <= 0) && verbose
	if(is.na(.test))
	cat("\nNote: missing values in density\n")	#
	else if(.test)
	cat("\nNote: negative values in density\n")	#
	#
	if(no.censor){
		o.B <- numeric(0)
	}else{
		#
		## compute for censorings
		mu0.link <- link.mu(covar.mu0 %*% x[2:(1 + ncovar.mu)])
		c0.link <- link.c(covar.c0 %*% x[(2 + ncovar.mu):(1 + ncovar.mu + ncovar.c)])
		if(protect)
			c0.link <- pmax(c0.link, 0)
		if(!use.grad){
			o.B <- f.B(t = .time0, mu = mu0.link, c1 = c0.link, tau = x[1])
		}else{
			o.B <- f.B.neglog.deriv(t = .time0, mu = mu0.link, c1 = c0.link, tau = x[1])
		}
		#
		.test <- any(o.B <= 0) && verbose
		if(is.na(.test))
			cat("\nNote: missing values in survival function\n")
		else if(.test)
			cat("\nNote: negative values in survival function\n")
	}
	#
	i <<- i+1 # ASSIGN TO ENCLOSURE
	if(verbose) cat(i, " ")
	return(list(o.fb = o.fb, o.B = o.B))	#
}# END .f.bB
##
## COMPUTE LIKELIHOOD
##
nloglike <- function(x){
	##
	## COMPUTE NEGATIVE LOG-LIKELIHOOD, WITHOUT DERIVATIVES
	.res <- .f.bB(x, use.grad = F)
	o.fb <- .res$o.fb
	o.B <- .res$o.B
	.ut <- c(o.fb, o.B)
	# IMPORTANT: WHEN USING GRADIENT, negative log IS ALREADY TAKEN!
	.ut <- - sum(log(.ut))
	return(.ut)	#
}# END nloglike
#
nloglike.grad <- function(x){
	##
	## COMPUTE NEGATIVE LOG-LIKELIHOOD, WITH DERIVATIVES
	.res <- .f.bB(x, use.grad = T)
	o.fb <- .res$o.fb
	o.B <- .res$o.B
	.ut <- c(o.fb, o.B)
	#
	# IMPORTANT: WHEN USING GRADIENT, negative log IS ALREADY TAKEN!
	.ut <- list(nll = sum(.ut))
	##
	## COMPUTE NEGATIVE LOG-LIKELIHOOD, WITH DERIVATIVES
	#
	gradient.fb <- attr(o.fb, "gradient")
	gradient.B <- attr(o.B, "gradient")	#
	gradient <- rbind(gradient.fb, gradient.B)
	#
	mu10.link.deriv <- as.numeric(link.mu.deriv(covar.mu.10 %*% x[2:(1 + ncovar.mu)])) #
	c10.link.deriv <- as.numeric(link.c.deriv(covar.c.10 %*% x[(2 + ncovar.mu):(1 + ncovar.mu + ncovar.c)]))
	#
	gradient <- cbind(gradient[, "tau"], mu10.link.deriv * covar.mu.10 * gradient[, "mu"] , c10.link.deriv * covar.c.10 * gradient[, "c1"])	#
	gradient <- colSums(gradient)
	names(gradient) <- termnames
	#
if(use.hess){
#		# EXTEND HESSIAN TO ALL COVARIATE PARAMETERS, NOT ONLY c1
#		# ikke klargjort for hess over mu
#				hessian.fb <- attr(o.fb, "hessian")
#				hessian.B <- attr(o.B, "hessian")
#				# disse to har dim n1 x 3 x 3 og n0 x 3 x 3
#
#				hessian <- f.join(hessian.fb, hessian.B)	#
#				# denne har dim n1+n0 x 3 x 3
#
#		#		f.vis(str(hessian))
#		#		stop("xxx")
#				# denne maa settes opp ogsaa for mu:
#				hessian <- hessian[, c("mu", "tau", rep("c1", ncovar)), c("mu", 
#					"tau", rep("c1", ncovar))] * covar.10.hessian	#
#		#
		}
	.ut$gradient <- gradient
#		if(use.hess) .ut$hessian <- hessian	#
	return(.ut)
}# END nloglike.grad
#
# GET STARTING VALUES:
grovest <- f.grovest3(formula = formula.mu, data = data, mu = mu, tau = tau, init.plot = F)	#
startverdier <- c(tau = grovest["est.tau"], c(grovest["est.mu"], rep(0, ncovar.mu - 1)), grovest["est.c"], rep(0, ncovar.c - 1))	#
names(startverdier) <- termnames # doesn't work since (Intercept) is not a valid name
#
## MERK: STARTVERDIER KAN BEREGNES SMARTERE ETTER BRUK AV f.grovest3 (SE PAA GJENNOMSNITT AV KOVARIATER, LOG ETC....)
##
## COMPUTE ESTIMATES:
##
#
cat("Running optimization....\n")
#
if(!use.gradient) {
	res <- optimx(par = startverdier, fn = nloglike, method = opti.method, control = list(reltol = 1e-14, maxit = 100000), hessian = T)
}else{
	# note that using nloglike - not nloglike.grad - for computing fn, saves time
	.f.gradient <- function(x) nloglike.grad(x)$gradient
	res <- optimx(par = startverdier, fn = nloglike, gr = .f.gradient, method = opti.method, control = list(reltol = 1e-14, maxit = 100000), hessian = T)
}
#
if(res$convcode != 0) cat("Possible problem with convergence!\n")


#		if(F & use.gradient) {
#		## hvorfor trengtes dett?
#			nll.deriv <- function(x)
#			{
#				temp <- nloglike.grad(x)	#
#				gradient <- rep(1, dim(covar.c.10)[1]) %*% attr(
#				  temp, "gradient")
#				hessian <- apply(attr(temp, "hessian"), c(2, 3),
#				  sum)
#				list(gradient = gradient, hessian = hessian[row(
#				  hessian) <= col(hessian)])
#			}
#			res <- nlminb(objective = nll, gradient = nll.deriv, 
#				hessian = T, start = startverdier)
#		}
#
## EXTRACT
.npar <- attr(res, "npar")
if(.npar != 1 + ncovar.mu + ncovar.c) stop("Something strange with the number of parameters!") ## IKKE EGENTLIG NODVENDIG
parameters <- unlist(res[1, 1:.npar]) # need unlist since data frame
# names(parameters) <- termnames
parameters["tau"] <- abs(parameters["tau"])	# unconstrained in estimation, should be positive
vcovmat <- solve(attr(res, "details")[1,]$nhatend)
dimnames(vcovmat) <- list(termnames, termnames)
max.loglike <-  - res$value
.AIC <- - 2 * max.loglike + 2 * .npar
#


## hva trengs denne for? ref. til grovest??
centered.intercept <- parameters["Intercept.c"]	#
#

#
out <- list(coefficients = parameters, cov.unscaled = vcovmat, loglik = max.loglike, AIC = .AIC, call = match.call(), link.mu = link.mu, link.c = link.c, opti.method = opti.method, centered.intercept = centered.intercept, result = res)
class(out) <- "invGauss"
return(out)
}

