#' Estimation of polygenic model
#' 
#' Estimation of polygenic model using a hierarchical generalized linear model 
#' (HGLM; Lee and Nelder 1996. \code{hglm} package; Ronnegard et al. 2010).
#' 
#' The algorithm gives extended quasi-likelihood (EQL) estimates
#' for restricted maximum likelihood (REML) (Ronnegard et al. 2010).
#' Implemented is the inter-connected GLM interpretation of HGLM (Lee and Nelder 2001).
#' Testing on fixed and random effects estimates are directly done using inter-connected
#' GLMs. Testing on dispersion parameters (variance components) can be done by extracting
#' the profile likelihood function value of REML.
#' 
#' 
#' @param formula Formula describing fixed effects to be used in analysis, e.g. 
#' y ~ a + b means that outcome (y) depends on two covariates, a and b. 
#' If no covariates used in analysis, skip the right-hand side of the 
#' equation.
#' @param kinship.matrix Kinship matrix, as provided by e.g. ibs(,weight="freq"), 
#' or estimated outside of GenABEL from pedigree data.
#' @param data An (optional) object of \code{\link{gwaa.data-class}} or a data frame with 
#' outcome and covariates.
#' @param family a description of the error distribution and link function to be 
#' used in the mean part of the model (see \code{\link{family}} for details of 
#' family functions).
#' @param conv 'conv' parameter of \code{hglm} (stricter than default, for great precision, use 1e-8).
#' @param maxit 'maxit' parameter of \code{hglm} (stricter than default).
#' @param ... other parameters passed to \code{hglm} call
#' 
#' @author Xia Shen, Yurii Aulchenko
#' 
#' @references 
#' Ronnegard, L, Shen, X and Alam, M (2010). hglm: A Package For Fitting 
#' Hierarchical Generalized Linear Models. \emph{The R Journal}, \bold{2}(2), 20-28.
#'
#' Lee, Y and Nelder JA (2001). Hierarchical generalised linear models: 
#' A synthesis of generalised linear models, random-effect models 
#' and structured dispersions. \emph{Biometrika}, \bold{88}(4), 987-1006.
#'
#' Lee, Y and Nelder JA (1996). Hierarchical Generalized Linear Models. 
#' \emph{Journal of the Royal Statistical Society. Series B}, \bold{58}(4), 619-678.
#' 
#' @seealso 
#' \code{\link{polygenic}},
#' \code{\link{mmscore}},
#' \code{\link{grammar}}
#' 
#' @examples 
#' data(ge03d2ex.clean)
#' set.seed(1)
#' df <- ge03d2ex.clean[,sample(autosomal(ge03d2ex.clean),1000)]
#' gkin <- ibs(df,w="freq")
#' 
#' # ----- for quantitative traits
#' h2ht <- polygenic_hglm(height ~ sex + age, kin = gkin, df)
#' # ----- estimate of heritability
#' h2ht$esth2
#' # ----- other parameters
#' h2ht$h2an
#' 
#' # ----- test the significance of fixed effects 
#' # ----- (e.g. quantitative trait)
#' h2ht <- polygenic_hglm(height ~ sex + age, kin = gkin, df)
#' # ----- summary with standard errors and p-values
#' summary(h2ht$hglm)
#' # ----- output for the fixed effects part
#' # ...
#' # MEAN MODEL
#' # Summary of the fixed effects estimates 
#' #              Estimate Std. Error t value Pr(>|t|)    
#' # (Intercept) 169.53525    2.57624  65.807  < 2e-16 ***
#' # sex          14.10694    1.41163   9.993  4.8e-15 ***
#' # age          -0.15332    0.05208  -2.944  0.00441 ** 
#' # ---
#' # Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#' # Note: P-values are based on 69 degrees of freedom
#' # ...
#' # ----- extract fixed effects estimates and standard errors
#' h2ht$h2an
#' 
#' # ----- test the significance of polygenic effect
#' nullht <- lm(height ~ sex + age, df)
#' l1 <- h2ht$ProfLogLik
#' l0 <- as.numeric(logLik(nullht))
#' # the likelihood ratio test (LRT) statistic
#' (S <- -2*(l0 - l1))
#' # 5% right-tailed significance cutoff 
#' # for 50:50 mixture distribution between point mass 0 and \chi^{2}(1)
#' # Self, SG and Liang KY (1987) Journal of the American Statistical Association.
#' qchisq(((1 - .05) - .50)/.50, 1)
#' # p-value
#' pchisq(S, 1, lower.tail = FALSE)/2
#' 
#' # ----- for binary traits
#' h2dm <- polygenic_hglm(dm2 ~ sex + age, kin = gkin, df, family = binomial(link = 'logit'))
#' # ----- estimated parameters
#' h2dm$h2an
#' 
#' 
#' @keywords htest
#' 
"polygenic_hglm" <- function(formula, kinship.matrix, data, family = gaussian(), conv = 1e-6, maxit = 100, ...)
{
	if (!require(hglm))
		stop("this function requires 'hglm' package to be installed")
	if (!missing(data)) if (is(data,"gwaa.data")) 
		{
			checkphengen(data)
			data <- phdata(data)
		}
	if (!missing(data)) 
		if (!is(data,"data.frame")) 
			stop("data should be of gwaa.data or data.frame class")
	allids <- data$id
	
	relmat <- kinship.matrix
	relmat[upper.tri(relmat)] <- t(relmat)[upper.tri(relmat)]
	
	call <- match.call()
	if (inherits(try(length(formula), silent = TRUE), "try-error")) 
		formula <- as.formula(paste(as.character(call)[2], '~ 1'))
	
	mf <- model.frame(formula, data, na.action = na.omit, drop.unused.levels = TRUE)
	y <- model.response(mf)
	desmat <- model.matrix(formula, mf)
	phids <- rownames(data)[rownames(data) %in% rownames(mf)]
	mids <- (allids %in% phids)
	relmat <- relmat[mids,mids]
	relmat <- relmat*2.0
	s <- svd(relmat)
	L <- s$u %*% diag(sqrt(s$d))
	#res_hglm <- try( hglm(y = y, X = desmat, Z = L, family = family, conv = conv, maxit = maxit, ... ) )
	res_hglm <- hglm(y = y, X = desmat, Z = L, family = family, conv = conv, maxit = maxit, ... )
	#sum_res_hglm <- summary(res_hglm)
	
	#if (is(try(tmp <- res_hglm$fixef),"try-error")) {
	#	warning("HGLM failed")
	#	return(res_hglm)
	#}
	#if (length(names(res_hglm$fixef))<1) 
	#	names(res_hglm$fixef) <- paste("fe",1:length(res_hglm$fixef),sep="")
	
	out <- list()
	
	out$measuredIDs <- mids
	out$hglm <- res_hglm
	out$h2an <- list()
	tVar <- res_hglm$varRan + res_hglm$varFix
	out$esth2 <- res_hglm$varRan/tVar
	out$h2an$estimate <- c(res_hglm$fixef,out$esth2, tVar)
	names(out$h2an$estimate)[(length(out$h2an$estimate) - 1):(length(out$h2an$estimate))] <- 
			c("h2","totalVar")
	if (is.null(res_hglm$SeFe)) {
		warning("Convergence failure HGLM!")
		out$h2an$se <- rep(NA,length(res_hglm$fixef)+2)
	} else {
		out$h2an$se <- c(res_hglm$SeFe, NA, NA)
	}
	names(out$h2an$se) <- 
			c(names(res_hglm$fixef), "h2","totalVar")
	out$pgresidualY <- rep(NA, length(mids))
	out$pgresidualY[mids] <- y - res_hglm$fv
	out$residualY <- rep(NA, length(mids))
	out$residualY[mids] <- y - desmat %*% res_hglm$fixef
	Sigma <- tVar*out$esth2*relmat + diag(tVar*(1 - out$esth2), ncol = length(y), nrow = length(y))
	out$InvSigma <- ginv(Sigma)
	#out$InvSigma <- ginv(as.matrix(res_hglm$vcov))
	#out$ProfLogLik <- res_hglm$ProfLogLik
	
	logdetSigma <- sum(log(eigen(Sigma, only.values = TRUE)$values))
	temp <- determinant(t(desmat)%*%out$InvSigma%*%desmat, logarithm = TRUE)
	out$ProfLogLik <- as.numeric(- .5*logdetSigma 
					- .5*t(out$residualY[mids])%*%out$InvSigma%*%out$residualY[mids]
					- .5*temp$modulus*temp$sign)
	
	class(out) <- "polygenic"
	out
}
