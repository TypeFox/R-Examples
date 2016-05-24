BayesSAE <-
function(formula, innov = "normal", df = NULL, b = NA, spatial = FALSE, tran = "F", prox = 
     NULL, beta.start = NULL, theta.start = NULL, lam.start = runif(1), prior = NULL, mcmc = 
     5000, burnin = 2500, thin = 5, data){
     ## specify the model
     call <- match.call()
     if (missing(data))
         data <- environment(formula)
   	 mf <- match.call(expand.dots = FALSE)
   	 m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
   	 mf <- mf[c(1L, m)]
   	 mf$drop.unused.levels <- TRUE
   	 oformula <- as.formula(formula)
   	 formula <- as.Formula(formula)
   	 if (length(formula)[2L] < 2L){
         formula <- as.Formula(formula(formula), ~1)	
     }	
     else {
         if (length(formula)[2L] > 2L) {
             formula <- Formula(formula(formula, rhs = 1:2))
             warning("formula must not have more than two RHS parts")
         }
     }
     mf$formula <- formula
     mf[[1L]] <- as.name("model.frame")
     mf <- eval(mf, parent.frame())
     mt <- terms(formula, data = data)
     mtX <- terms(formula, data = data, rhs = 1L)
     mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
     Y <- model.response(mf, "numeric")
     X <- model.matrix(mtX, mf)
     Z <- model.matrix(mtZ, mf)
     Z <- as.vector(Z[,2])
     m <- length(Y)
     p <- ncol(X)
	 
     ## identify the weights
     if (is.na(b))
         b <- rep(1, m)

     ## determine the initial values	 
     if (tran != "F" && tran != "log" && tran != "logit")
         stop("Only log or logit transformations are allowed")
     if (is.null(theta.start) || is.null(beta.start)){
         if (tran == "F"){
             if (is.null(beta.start))
	             beta.start = as.vector(coef(lm(Y~X[,-1], weights = 1 / b)))
             if (is.null(theta.start))
                 theta.start = as.vector(Y)			 
         }
	     else if (tran == "log"){
             if (is.null(beta.start))
                 beta.start = as.vector(coef(lm(log(Y)~X[,-1], weights = 1 / b)))
             if (is.null(theta.start))
                 theta.start = as.vector(Y)
         }
         else {
             if (is.null(beta.start))
                 beta.start = as.vector(coef(lm(log(Y)-log(1-Y)~X[,-1], weights = 1 / b)))
             if (is.null(theta.start))
                 theta.start = as.vector(Y)
         }		 
     }

     ## check the initial values
     if (tran == "log" && any(theta.start <= 0))
         stop("initial values for theta's must be positive")
     if (tran == "logit" && (any(theta.start <= 0) || any(theta.start >= 1)))
         stop("initial values for theta's must be bounded between 0 and 1")
     if (spatial){ 
         if ((lam.start <= 0)||(lam.start >= 1))
             stop("initial value for lambda must be bounded between 0 and 1")
     }
	 
     ## check the response variable
     if (tran == "log" && any(Y <= 0))
         stop("response variable must be positive")
     if (tran == "logit" && (any(Y <= 0) || any(Y >= 1)))
         stop("response variable must bounded between 0 and 1")	 
     
     ## check the distribution of innovation
     if (innov != "normal" && innov != "t")
         stop("innovations in the linking model must be normally or t distributed")
 
     ## check the spatial model
     if (spatial && tran != "F")
         stop("unmatched spatial models are not supported currently")
		 
     ## obtain the prior
     if (is.null(prior)){
         beta.type <- "non_in"
         beta.prior <- 0
         sigv.type <- "unif"
         sigv.prior <- list(eps2 = 0.0001)
     }
     else{
         beta.type <- prior$beta.type
         beta.prior <- prior$beta.prior
         sigv.type <- prior$sigv.type
         sigv.prior <- prior$sigv.prior	     
     }	 
	 if (innov == "t"){
         if (is.null(prior))
             sig2.prior <- list(ai = rep(0.05, m), bi = rep(0.05, m))
         else
             sig2.prior <- prior$sig2.prior
     }     
     ## check prior for beta
     if (beta.type != "normal" && beta.type != "non_in")
         stop("only normal or non-informative distribution is allowed as prior for beta's at present")
     if (sigv.type != "inv_gamma" && sigv.type != "unif")
         stop("only invert gamma or uniform distribution is allowed as prior for sigma2v at present")
     if (beta.type == "normal")
         beta.type <- 0
     else
         beta.type <- 1
     if (beta.type == 0){
         beta0 <- beta.prior$beta0
         if (length(beta0) != p)
             stop("number of regressors and length of beta0 are inconsistent")
         eps1 <- beta.prior$eps1
         if (eps1 <= 0)
             stop("eps1 must be positive")
         beta.prior <- c(beta0, eps1)
     }
     else
         beta.prior <- 0
     ## check prior for sigma2v
     if (sigv.type == "inv_gamma")
         sigv.type <- 0
     else
         sigv.type <- 1
     if (sigv.type == 0){
         a0 <- sigv.prior$a0
         b0 <- sigv.prior$b0
         if (a0 <= 0 || b0 <= 0)
             stop("a0 and b0 must be positive")
         sigv.prior <- c(a0 ,b0)
     }
     else{
         sigv.prior <- sigv.prior$eps2
         if (sigv.prior <= 0)
             stop("eps2 must be positive")
     }
     ## check prior for sigma2
     if (innov == "t"){
         ai <- sig2.prior$ai
         bi <- sig2.prior$bi
         if (any(ai <= 0) || any(bi <= 0))
             stop("elements in ai and bi must be positive")
         if (length(ai) != m || length(bi) != m)
             stop("lengths of ai and bi should be the same as that of domains")
     }
	 
     ## check values of degree of freedom
     if (innov == "t"){
         if (any(df <= 0))
             stop("elements in df must be positive")
         if (length(df) != m)
             stop("length of df should be the same as that of domains")
     }     
     
     ## calculate subset
     subset <- seq(from = burnin+1, to = mcmc, by = thin)	 
	 
     ## begin sampling	 
     if (spatial){
         li1 <- prox[,1] 
         li2 <- prox[,2] 
         num <- rep(0, m)
         for (i in 1:m){
             num[i] = sum(li1 == i) + sum(li2 == i)
         }
         if (innov == "normal")
             result <- BayesSFH(theta.start, beta.start, lam.start, Y, t(X), Z, li1, li2, 
                 num, mcmc, beta.prior, sigv.prior, beta.type, sigv.type)
         else
             result <- BayesSYC(theta.start, beta.start, lam.start, Y, t(X), Z, li1, li2, 
                 num, mcmc, beta.prior, sigv.prior, c(ai, bi, df), beta.type, sigv.type)
     }
     else{
         if (innov == "normal"){
             if (tran == "F")
                 result <- BayesFH(theta.start, beta.start, Y, t(X), b, Z, mcmc, 
                     beta.prior, sigv.prior, beta.type, sigv.type)
             else if (tran == "log")
                 result <- BayesUFH(theta.start, beta.start, Y, t(X), b, Z, mcmc, 
                     beta.prior, sigv.prior, beta.type, 1, sigv.type)
             else
                 result <- BayesUFH(theta.start, beta.start, Y, t(X), b, Z, mcmc, 
                     beta.prior, sigv.prior, beta.type, 2, sigv.type)
         }
         else{
             if (tran == "F")
                 result <- BayesYC(theta.start, beta.start, Y, t(X), b, Z, mcmc, 
                     beta.prior, sigv.prior, c(ai, bi, df), beta.type, sigv.type)
             else if (tran == "log")
                 result <- BayesUYC(theta.start, beta.start, Y, t(X), b, Z, mcmc, 
                     beta.prior, sigv.prior, c(ai, bi, df), beta.type, 1, sigv.type)
             else
                 result <- BayesUYC(theta.start, beta.start, Y, t(X), b, Z, mcmc, 
                     beta.prior, sigv.prior, c(ai, bi, df), beta.type, 2, sigv.type)
         }
     }
     result$innov <- innov
     result$prox <- prox
     result$X <- X
     result$Y <- Y
     result$Z <- Z
     result$b <- b
     result$m <- m
     result$p <- p
     result$subset <- subset
     theta <- t(result$theta)
     beta <- t(result$beta)
     sigv <- result$sigv
     if (!(is.null(result$lam)))
         lam <- result$lam
     if (!(is.null(result$sig2)))
         sig2 <- t(result$sig2)
     if (spatial){
         li1 <- prox[,1] 
         li2 <- prox[,2] 
         num <- rep(0, m)
         for (i in 1:m){
             num[i] = sum(li1 == i) + sum(li2 == i)
         }
         if (innov == "normal"){
             result$type <- "SFH"
             HB <- theta.HB(result, subset)
             criter <- model.compare(result, subset, poest = "mean")
             result <- list(mcmc = mcmc(data.frame(theta = theta, beta = beta, 
                 sigv = sigv, lam = lam), start = burnin+1, end = mcmc, thin = thin), 
                 lambda.rate = result$lam.rate, D_avg = criter$D_avg, D_theta.hat = 
                 criter$D_theta.hat, DIC = criter$DIC, theta.HB = HB)
             result$type <- "SFH"
             result$HB <- HB
         }
         else{
             result$type <- "SYC"
             HB <- theta.HB(result, subset)
             criter <- model.compare(result, subset, poest = "mean")
             result <- list(mcmc = mcmc(data.frame(theta = theta, beta = beta, 
                 sigv = sigv, sig2 = sig2, lam = lam), start = burnin+1, end = mcmc, 
                 thin = thin), lam.rate = result$lam.rate, D_avg = criter$D_avg, D_theta.hat
                 = criter$D_theta.hat, DIC = criter$DIC, theta.HB = HB)	
             result$type <- "SYC"
             result$HB <- HB
         }				 
     }
     else{
         if (innov == "normal"){
             if (tran == "F"){
                 result$type <- "FH"
                 HB <- theta.HB(result, subset)
                 criter <- model.compare(result, subset, poest = "mean")
                 result <- list(mcmc = mcmc(data.frame(theta = theta, beta = beta, 
                     sigv = sigv), start = burnin+1, end = mcmc, thin = thin), D_avg = 
                     criter$D_avg, D_theta.hat = criter$D_theta.hat, DIC = criter$DIC, 
                     theta.HB = HB)
                 result$type <- "FH"
                 result$HB <- HB
             }
             else if (tran == "log"){
                 criter <- model.compare(result, subset, poest = "mean")
                 result <- list(mcmc = mcmc(data.frame(theta = theta, beta = beta, 
                     sigv = sigv), start = burnin+1, end = mcmc, thin = thin), theta.rate
                     = result$theta.rate, D_avg = criter$D_avg, D_theta.hat = criter$D_theta.hat, 
                     DIC = criter$DIC)
                 result$type <- "UFH"
             }
             else{
                 criter <- model.compare(result, subset, poest = "mean")
                 result <- list(mcmc = mcmc(data.frame(theta = theta, beta = beta, 
                     sigv = sigv), start = burnin+1, end = mcmc, thin = thin), theta.rate = 
                     result$theta.rate, D_avg = criter$D_avg, D_theta.hat = criter$D_theta.hat, 
                     DIC = criter$DIC)
                 result$type <- "UFH"
             }
         }
         else{
             if (tran == "F"){
                 result$type <- "YC"
                 HB <- theta.HB(result, subset)
                 criter <- model.compare(result, subset, poest = "mean")
                 result <- list(mcmc = mcmc(data.frame(theta = theta, beta = beta, 
                     sigv = sigv, sig2 = sig2), start = burnin+1, end = mcmc, thin = thin), 
                     D_avg = criter$D_avg, D_theta.hat = criter$D_theta.hat, DIC = criter$DIC, 
                     theta.HB = HB)
                 result$type <- "YC"
                 result$HB <- HB
             }
             else if (tran == "log"){
                 criter <- model.compare(result, subset, poest = "mean")
                 result <- list(mcmc = mcmc(data.frame(theta = theta, beta = beta, 
                     sigv = sigv, sig2 = sig2), start = burnin+1, end = mcmc, thin = thin), 
					 theta.rate = result$theta.rate, D_avg = criter$D_avg, 
					 D_theta.hat = criter$D_theta.hat, DIC = criter$DIC)
                 result$type <- "UYC"					 
             }
             else{
                 criter <- model.compare(result, subset, poest = "mean")
                 result <- list(mcmc = mcmc(data.frame(theta = theta, beta = beta, 
                     sigv = sigv, sig2 = sig2), start = burnin+1, end = mcmc, thin = thin), 
                     theta.rate = result$theta.rate, D_avg = criter$D_avg, D_theta.hat = 
                     criter$D_theta.hat, DIC = criter$DIC)
                 result$type <- "UYC"
             }
         }
     }		 
     result$innov <- innov
     result$prox <- prox
     result$X <- X
     result$Y <- Y
     result$Z <- Z
     result$b <- b
     result$m <- m
     result$p <- p
     result$mf <- mf
     result$spatial <- spatial
     result$tran <- tran	 
     result$subset <- subset
     result$call <- call
     class(result) <- "BayesSAE"
     result
}
