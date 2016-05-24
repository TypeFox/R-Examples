#' Proportional Hazards Model with Mixed Effects
#' 
#' Fits a proportional hazards regression model incorporating random effects.
#' The function implements an EM algorithm using Markov Chain Monte Carlo
#' (MCMC) at the E-step as described in Vaida and Xu (2000).
#' 
#' The proportional hazards model with mixed effects is equipped to handle
#' clustered survival data. The model generalizes the usual frailty model by
#' allowing log-linearl multivariate random effects. The software can only
#' handle random effects from a multivariate normal distribution. Maximum
#' likelihood estimates of the regression parameters and variance components is
#' gotten by EM algorithm, with Markov chain Monte Carlo (MCMC) used in the
#' E-step.
#' 
#' Care must be taken to ensure the MCMC-EM algorithm has converged, as the
#' algorithm stops after MAXSTEP iterations. No convergence criteria is
#' implemented. It is advised to plot the estimates at each iteration using the
#' \code{\link[phmm]{plot.phmm}} method. For more on MCMC-EM convergence see
#' Booth and Hobart (1999).
#' 
#' @param formula model formula for the fixed and random components of the
#' model (as in \code{\link[lme4]{lmer}}). An intercept is implicitly included
#' in the model by estimation of the error distribution. As a consequence
#' \code{-1} in the model formula does not have any effect.  The left-hand side
#' of the \code{formula} must be a \code{\link[survival]{Surv}} object.
#' @param data optional data frame in which to interpret the variables occuring
#' in the formulas.
#' @param subset subset of the observations to be used in the fit.
#' @param na.action function to be used to handle any \code{NA}s in the data.
#' The user is discouraged to change a default value \code{na.fail}.
#' @param Sigma initial covariance matrix for the random effects. Defaults to
#' "identity".
#' @param varcov constraint on \code{Sigma}. Currently only \code{"diagonal"}
#' is supported.
#' @param NINIT number of starting values supplied to Adaptive Rejection
#' Metropolis Sampling (ARMS) algorithm.
#' @param VARSTART starting value of the variances of the random effects.
#' @param MAXSTEP number of EM iterations.
#' @param CONVERG iteration after which Gibbs sampling size changes from Gbs to
#' Gbsvar.
#' @param Gbs initial Gibbs sampling size (until CONVERG iterations).
#' @param Gbsvar Gibbs sampling size after CONVERG iterations.
#' @param verbose Set to \code{TRUE} to print EM steps.
#' @param maxtime maximum time in seconds, before aborting EM iterations.
#' Defaults to 120 seconds.
#' @param random The argument \code{random} is no longer used. Random
#' components are are expressed in \code{formula}.
#' @return The function produces an object of class "phmm" consisting of:
#' @returnItem steps a matrix of estimates at each EM step;
#' @returnItem bhat empirical Bayes estimates of expectation of random effects;
#' @returnItem sdbhat empirical Bayes estimates of standard deviation of random
#' effects;
#' @returnItem coef the final parameter estimates for the fixed effects;
#' @returnItem var the estimated variance-covariance matrix;
#' @returnItem loglik a vector of length four with the conditional
#' log-likelihood and marginal log-likelihood estimated by Laplace
#' approximation, reciprocal importance sampling, and bridge sampling (only
#' implemented for \code{nreff} < 3);
#' @returnItem lambda the estimated baseline hazard;
#' @returnItem Lambda the estimated cumulative baseline hazard.
#' @seealso \code{\link[survival]{survfit}}, \code{\link[survival]{Surv}}.
#' @references Gilks WR and Wild P. (1992) Adaptive rejection sampling for
#' Gibbs sampling. Applied Statistics 41, pp 337-348.
#' 
#' Donohue, MC, Overholser, R, Xu, R, and Vaida, F (January 01, 2011).
#' Conditional Akaike information under generalized linear and proportional
#' hazards mixed models. \emph{Biometrika}, 98, 3, 685-700.
#' 
#' Vaida F and Xu R. 2000. "Proportional hazards model with random effects",
#' \emph{Statistics in Medicine,} 19:3309-3324.
#' 
#' Gamst A, Donohue M, and Xu R (2009). Asymptotic properties and empirical
#' evaluation of the NPMLE in the proportional hazards mixed-effects model.
#' Statistica Sinica, 19, 997.
#' 
#' Xu R, Gamst A, Donohue M, Vaida F, and Harrington DP. 2006. Using Profile
#' Likelihood for Semiparametric Model Selection with Application to
#' Proportional Hazards Mixed Models. \emph{Harvard University Biostatistics
#' Working Paper Series,} Working Paper 43.
#' 
#' Booth JG and Hobert JP. Maximizing generalized linear mixed model
#' likelihoods with an automated Monte Carlo EM algorithm. \emph{Journal of the
#' Royal Statistical Society}, Series B 1999; 61:265-285.
#' @keywords survival
#' @examples
#' 
#' n <- 50      # total sample size
#' nclust <- 5  # number of clusters
#' clusters <- rep(1:nclust,each=n/nclust)
#' beta0 <- c(1,2)
#' set.seed(13)
#' #generate phmm data set
#' Z <- cbind(Z1=sample(0:1,n,replace=TRUE),
#'            Z2=sample(0:1,n,replace=TRUE),
#'            Z3=sample(0:1,n,replace=TRUE))
#' b <- cbind(rep(rnorm(nclust),each=n/nclust),rep(rnorm(nclust),each=n/nclust))
#' Wb <- matrix(0,n,2)
#' for( j in 1:2) Wb[,j] <- Z[,j]*b[,j]
#' Wb <- apply(Wb,1,sum)
#' T <- -log(runif(n,0,1))*exp(-Z[,c('Z1','Z2')]%*%beta0-Wb)
#' C <- runif(n,0,1)
#' time <- ifelse(T<C,T,C)
#' event <- ifelse(T<=C,1,0)
#' mean(event)
#' phmmd <- data.frame(Z)
#' phmmd$cluster <- clusters
#' phmmd$time <- time
#' phmmd$event <- event
#' 
#' fit.phmm <- phmm(Surv(time, event) ~ Z1 + Z2 + (-1 + Z1 + Z2 | cluster), 
#'    phmmd, Gbs = 100, Gbsvar = 1000, VARSTART = 1,
#'    NINIT = 10, MAXSTEP = 100, CONVERG=90)
#' summary(fit.phmm)
#' plot(fit.phmm)
#' 
#' @export phmm
phmm <- function (formula, data, subset, 
	na.action = na.fail, Sigma = "identity", varcov = "diagonal", 
	NINIT = 10, VARSTART = 1, MAXSTEP = 100, CONVERG = 90, Gbs = 100, 
	Gbsvar = 1000, verbose = FALSE, maxtime = 120, random)
{
    Call <- match.call()

    if (!missing(random)) {
        warning("The 'random' argument of phmm is deprecated")
        if (class(random) != 'formula' || length(random) !=2) 
            stop("Invalid random formula")
        j <- length(formula)   #will be 2 or 3, depending on if there is a y

        # Add parens to the random formula and paste it on
        formula[[j]] <- call('+', formula[[j]], call('(', random[[2]]))  
        }

    temp <- call('model.frame', formula= subbar(formula))
    for (i in c('data', 'subset', 'na.action'))
        if (!is.null(Call[[i]])) temp[[i]] <- Call[[i]]
    if (is.R()) m <- eval.parent(temp)
    else        m <- eval(temp, sys.parent())
        Y <- model.extract(m, "response")
        n <- nrow(Y)
        if (!inherits(Y, "Surv")) stop("Response must be a survival object")
        type <- attr(Y, "type")
        if (type!='right')
            stop(paste("phmm doesn't support '", type,
                              "' survival data", sep=''))

        # Check for penalized terms; the most likely is pspline
        pterms <- sapply(m, inherits, 'coxph.penalty')
        if (any(pterms)) {
            stop("Penalized terms are not supported by phmm")
            }

        flist <- formula1(formula)
        if (hasAbar(flist$fixed))
            stop("Invalid formula: a '|' outside of a valid random effects term")

        special <- c("strata", "cluster")
        Terms <- terms(flist$fixed, special)
        attr(Terms,"intercept")<- 1  #Cox model always has \Lambda_0
        strats <- attr(Terms, "specials")$strata
        cluster<- attr(Terms, "specials")$cluster
        if (length(cluster)) {
            stop ("A cluster() statement is invalid in phmm")
            }
        if (length(strats)) {
            stop ("A strata() statement is not supported in phmm")
            }
        else Z <- model.matrix(Terms, m)[,-1,drop=F]
    
    
    random <- unlist(strsplit(as.character(flist$random[[1]])[2], "|", fixed=TRUE))
    W <- model.matrix(as.formula(paste("~", random[1])), m)[,,drop=F]
    cluster <- model.matrix(as.formula(paste("~", random[2])), m)[,-1,drop=F]
    if(ncol(cluster) != 1) stop("Unsupported cluster structure")
    nclust <- length(unique(cluster))
    n <- nrow(Z)
    nfixed <- ncol(Z)
	names.random <- colnames(W)

    nrandom <- ncol(W)
    if (nrandom ==0) stop("No random effects terms found")

	if(Sigma == "identity") Sigma0 = diag(1, nrandom)
	invSigma = rbind(0, cbind(0, solve(Sigma0)))
	detSigma = det(Sigma0)
	Sigma    = rbind(0, cbind(0, solve(Sigma0)))

	if(varcov == "diagonal"){ varcov = 2
	}else{ stop(paste("\nUnknown structure specified for var-covariance matrix:", 
		varcov))}

    if(verbose){
		cat("\nProportional Hazards Mixed-Effects Model fit with MCMC-EM\n")
		date0 <- date()
	}

	fit <- .C("phmm", 
		X = as.numeric(c(0,Y[,1])), 
		Z = rbind(0, as.matrix(Z)), 
		W = rbind(0, as.matrix(W)),
		delta = as.integer(c(0,Y[,2])), 
		cluster = as.integer(c(0,cluster)),
		varcov = as.integer(varcov),
		Sigma = Sigma,
		invSigma = invSigma,
		detSigma = as.double(detSigma),
		NINIT=as.integer(10),
		MAXSTEP=as.integer(MAXSTEP),
		maxtime=as.double(maxtime),
		CONVERG=as.integer(CONVERG),
		emstep = as.integer(0),
		Gbs = as.integer(Gbs),
		Gbsvar = as.integer(Gbsvar),
		n = as.integer(n),
		nfixed = as.integer(nfixed),
		nrandom = as.integer(nrandom),
		bhat = double(nclust*nrandom),
		sdbhat = double(nclust*nrandom),
		steps = double((MAXSTEP+1)*(nfixed+nrandom)),
		var = double((nfixed+nrandom+n+1)^2),
		verbose = as.integer(verbose),
		bridgeC = double(1),
		laplacetot = double(1),
		llaplace = double(1), 
		limport = double(1), 
		lbridge = double(1),
		lambda = double(n+1),
		Lambda = double(n+1),
		PACKAGE="phmm" )

	if(varcov == 2) fit$varcov = "diagonal"
	
	class(fit) <- "phmm"

	fit$Call <- Call
	fit$formula <- formula
	fit$random <- random
	fit$verbose <- as.logical(fit$verbose)
	fit$steps <- matrix(fit$steps, nrow = MAXSTEP+1, byrow = TRUE)
	colnames(fit$steps) <- c(colnames(Z), paste("var(", names.random, ")", sep = ""))
	rownames(fit$steps) <- 0:MAXSTEP
	
	fit$var <- matrix(fit$var, nrow = nfixed+nrandom+n+1, byrow = FALSE)
	fit$var <- fit$var[1 + 1:(nfixed+nrandom), 1 + 1:(nfixed+nrandom)]
	rownames(fit$var) <- colnames(fit$var) <- c(colnames(Z), names.random)
	fit$varFix <- fit$var[1:nfixed, 1:nfixed]
	
	fit$Sigma0 <- Sigma0
	fit$Sigma <- matrix(fit$Sigma, nrow = nrandom+1, byrow = FALSE)
	fit$Sigma <- fit$Sigma[1 + 1:nrandom, 1 + 1:nrandom]
	if(nrandom == 1){ names(fit$Sigma) = paste("var(", names.random, ")", sep = "")
	}else{
		colnames(fit$Sigma) = names.random
		rownames(fit$Sigma) = names.random
	}
	fit$invSigma <- matrix(fit$invSigma, nrow = nrandom+1, byrow = FALSE)
	fit$invSigma <- fit$invSigma[1 + 1:nrandom, 1 + 1:nrandom]

	fit$coefficients = fit$steps[MAXSTEP+1, 1:fit$nfixed]
	names(fit$coefficients) <- colnames(Z)
	fit$bhat <- matrix(fit$bhat, nrow = nclust, byrow = TRUE)
	fit$sdbhat <- matrix(fit$sdbhat, nrow = nclust, byrow = TRUE)
	
	rownames(fit$bhat) <- rownames(fit$sdbhat) <- unique(cluster) 
	colnames(fit$bhat) <- colnames(fit$sdbhat) <- colnames(W)	
	
	fit$X = Y[, 1]
	fit$Z = as.matrix(Z)
	fit$W = as.matrix(W)
	fit$delta = Y[, 2] 
	fit$cluster = cluster
	
	fit$Lambda = fit$Lambda[-1]
	fit$lambda = fit$lambda[-1]
	
	fit$bhat.long <- as.matrix(merge(
	   cbind(cluster = as.numeric(rownames(fit$bhat)), as.matrix(fit$bhat)), 
	   cbind(cluster = as.numeric(as.character(fit$cluster)))))[, -1]
	fit$Y <- Y
	
	if(fit$nrandom > 2) fit$lbridge <- NULL
	fit$loglik <- c(Conditional = loglik.cond(fit), 
		Laplace = fit$llaplace, RIS = fit$limport, BS = fit$lbridge)
		
	fit$linear.predictors <- linear.predictors(fit)
	fit$flist <- flist
	fit$Terms <- Terms

	fit <- fit[c('Call', 'Y', 'Z', 'W', 'cluster', 'varcov', 'Sigma', 
	   'invSigma', 'detSigma', 'NINIT', 'MAXSTEP', 'CONVERG', 'emstep', 
	   'Gbs', 'Gbsvar', 'n', 'nfixed', 'nrandom', 'bhat', 'bhat.long', 
	   'sdbhat', 'steps', 'var', 'verbose', 'loglik', 'lambda', 'Lambda', 
	   'varFix', 'coefficients', 'linear.predictors',
	   'formula', 'flist', 'Terms')]
	
	class(fit) <- "phmm"
    return(fit)
}

#' PHMM conditional log-likelihood
#' 
#' Function for computing log-likelihood conditional on the estimated random
#' effects from an object of class \code{phmm} returned by \code{phmm}.
#' 
#' 
#' @param x an object of class \code{phmm}.
#' @return The PHMM log-likelihood conditional on the estimated random effects.
#' @seealso \code{\link{phmm}}, \code{\link{phmm.cond.loglik}}
#' @keywords survival
loglik.cond <- function (x) UseMethod("loglik.cond")
loglik.cond.phmm <- function(x){
	#Function to compute conditional log-likelihood
	phmm.cond.loglik(time = x$Y[, 1], delta = x$Y[, 2], z = x$Z, beta = x$coef, w = x$W, b = as.matrix(x$bhat.long))
}

#' PHMM conditional log-likelihood
#' 
#' Function for computing log-likelihood conditional on the estimated random
#' effects from the data and specified parameter estimates of a PHMM.
#' 
#' 
#' @param time Follow-up time (right censored data).
#' @param delta The status indicator (0=alive, 1=dead; or \code{TRUE}=dead,
#' \code{FALSE}=alive).
#' @param z Numeric matrix (\code{N}x\code{nfixed}) of covariates for fixed
#' effects.
#' @param beta Fitted fixed effects coefficients (\code{p}-vector).
#' @param w Numeric matrix (\code{N}x\code{nrandom}) of covariates for random
#' effects.
#' @param b Numeric matrix (\code{N}x\code{nrandom}) of random effects
#' estimates.
#' @return The PHMM log-likelihood conditional on the estimated random effects.
#' @seealso \code{\link{phmm}}, \code{\link{loglik.cond}}
#' @keywords survival
phmm.cond.loglik <- function(time, delta, z, beta, w, b){
	#Function to compute conditional log-likelihood
    z <- as.matrix(z)
    wb <- matrix(0, nrow = nrow(w), ncol = ncol(w))
    for(i in 1:ncol(w)){
	  if(length(w[, i]) != length(b[, i])) stop("length(w[, i]) != length(b[, i])")
	  wb[, i] <- w[, i]*b[, i]
	  }
    if(!is.null(dim(wb))) wb <- apply(wb, 1, sum)
    numerator <- exp(z%*%beta+wb)
    denominator <- unlist(lapply(time, 
      FUN = function(x){
          sum(exp((z%*%beta+wb)[time >= x]))
          }))
    sum(ifelse(delta, 1, 0)*log(numerator/denominator))
}


#' Akaike Information Criterion for PHMM
#' 
#' Function calculating the Akaike information criterion for PHMM fitted model
#' objects, according to the formula \eqn{-2*log-likelihood +
#' k*rho}{-2*log-likelihood + k*npar}, where \eqn{npar}{npar} represents the
#' number of parameters in the fitted model. The function returns a list of AIC
#' calculations corresponding different likelihood estimations: conditional and
#' marginal likelihoods calculated by Laplace approximation, reciprocal
#' importance sampling, and bridge sampling (only implemented for nreff < 3).
#' The default k = 2, is for the usual AIC.
#' 
#' 
#' @aliases AIC.phmm
#' @param object a fitted PHMM model object of class \code{phmm},
#' @param ... optionally more fitted model objects.
#' @param k numeric, the penalty per parameter to be used; the default k = 2 is
#' the classical AIC.
#' @return Returns a list of AIC values corresonding to all available
#' log-likelihood values from the fit. See \code{\link{phmm}} for details of
#' the log-likelihood values.
#' @seealso \code{\link{phmm}}, \code{\link[stats]{AIC}}
#' @references Whitehead, J. (1980). Fitting Cox's Regression Model to Survival
#' Data using GLIM. Journal of the Royal Statistical Society. Series C, Applied
#' statistics, 29(3), 268-.
#' @keywords survival
#' @examples
#' 
#' n <- 50      # total sample size
#' nclust <- 5  # number of clusters
#' clusters <- rep(1:nclust,each=n/nclust)
#' beta0 <- c(1,2)
#' set.seed(13)
#' #generate phmm data set
#' Z <- cbind(Z1=sample(0:1,n,replace=TRUE),
#'            Z2=sample(0:1,n,replace=TRUE),
#'            Z3=sample(0:1,n,replace=TRUE))
#' b <- cbind(rep(rnorm(nclust),each=n/nclust),rep(rnorm(nclust),each=n/nclust))
#' Wb <- matrix(0,n,2)
#' for( j in 1:2) Wb[,j] <- Z[,j]*b[,j]
#' Wb <- apply(Wb,1,sum)
#' T <- -log(runif(n,0,1))*exp(-Z[,c('Z1','Z2')]%*%beta0-Wb)
#' C <- runif(n,0,1)
#' time <- ifelse(T<C,T,C)
#' event <- ifelse(T<=C,1,0)
#' mean(event)
#' phmmd <- data.frame(Z)
#' phmmd$cluster <- clusters
#' phmmd$time <- time
#' phmmd$event <- event
#' 
#' fit.phmm <- phmm(Surv(time, event) ~ Z1 + Z2 + (-1 + Z1 + Z2 | cluster), 
#'    phmmd, Gbs = 100, Gbsvar = 1000, VARSTART = 1,
#'    NINIT = 10, MAXSTEP = 100, CONVERG=90)
#' 
#' # Same data can be fit with lmer,
#' # though the correlation structures are different.
#' poisphmmd <- pseudoPoisPHMM(fit.phmm)
#' 
#' library(lme4)
#' fit.lmer <- lmer(m~-1+as.factor(time)+z1+z2+
#'   (-1+w1+w2|cluster)+offset(log(N)), 
#'   as.data.frame(as(poisphmmd, "matrix")), family=poisson)
#' 
#' fixef(fit.lmer)[c("z1","z2")]
#' fit.phmm$coef
#' 
#' VarCorr(fit.lmer)$cluster
#' fit.phmm$Sigma
#' 
#' logLik(fit.lmer)
#' fit.phmm$loglik
#' 
#' traceHat(fit.phmm)
#' 
#' summary(fit.lmer)@AICtab
#' AIC(fit.phmm)
#'
AIC.phmm <- function(object, ..., k = 2){
	if(object$varcov == "diagonal"){ 
		return(-2*object$loglik+k*(object$nrandom+object$nfixed))
	}else{ stop(paste("\nUnknown structure specified for var-covariance matrix:", 
		object$varcov))}
}
	
print.phmm <-
 function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nProportional Hazards Mixed-Effects Model fit by MCMC-EM\n")
  cat("  Model:", deparse(x$formula), "\n")
  cat("  Data:", deparse( x$Call$data ), "\n")
  cat("  Log-likelihood:\n")
	print(x$loglik[!is.null(x$loglik)], digits = digits, ...)
  cat("\nFixed effects:", deparse(x$flist$fixed), "\n")
  print(coef(x), digits = digits, ...)
  cat("\n")
  cat("Random effects:", deparse(x$flist$random[[1]]), "\n")
  cat("Estimated variance-covariance matrix:\n")
  print(x$Sigma, digits = digits, ...)
#  cat("Variance-Covariance:\n")
#	print(x$var, ...)
  cat("\nNumber of Observations:", x$n)
  cat("\nNumber of Groups: ", nrow(x$bhat))
  cat("\n\n")
}

print.summary.phmm <- print.phmm

summary.phmm <-
 function(object, ...)
{
	object$coefficients = cbind(Estimate = object$coef, 
	   Std.Error = sqrt(ifelse(diag(object$var)<0, NA, diag(object$var)))[1:object$nfixed])
	class(object)<-"summary.phmm"
	return(object)
}



#' Plots the convergence of MCMC-EM estimates from a PHMM
#' 
#' Plots the value of each parameter of the model at each iteration of the
#' MCMC-EM algorithm. For more on MCMC-EM convergence see Booth \& Hobart
#' (1999).
#' 
#' 
#' @param x \code{phmm} object return by \code{\link[phmm]{phmm}}
#' @param ... other arguments passed to \code{\link[lattice]{xyplot}}
#' @seealso \code{\link[phmm]{phmm}}.
#' @references Booth JG \& Hobert JP. Maximizing generalized linear mixed model
#' likelihoods with an automated Monte Carlo EM algorithm. \emph{Journal of the
#' Royal Statistical Society}, Series B 1999; 61:265-285.
#' @keywords survival
plot.phmm <-
 function(x, ...)
{
	x = as.data.frame(x$steps)
	colnames(x) = make.names(colnames(x))
	fm = paste(paste(colnames(x), collapse = ' + '), "EM.Step", sep = " ~ ")
	x$EM.Step = as.numeric(rownames(x))
	xyplot(formula(fm), data = x, type = "l", allow.multiple = TRUE, outer = TRUE, scales = "free", ...)
}
