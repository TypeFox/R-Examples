simglm2 <-
function(y, X, Z, link=1, beta, alpha, maxit, tol) {
##This program is used to produce the estimates of parameters for the case of independent 
##structure for Heterogeneity GLM method. model: y~Simplex(mu, Sigma), 
##logit(mu) = X * beta, log(Sigma_ij) = sigma0 + sigma1*z_1ij + sigma2*z_2ij

	if (link==1){
		gf <- function(x) gf1(x)
		hf <- function(x) hf1(x)
		gfd <- function(x) gfd1(x)
		}
	else if (link==2) {
		gf <- function(x) gf2(x)
		hf <- function(x) hf2(x)
		gfd <- function(x) gfd2(x)
		}
	else if (link==3) {
		gf <- function(x) gf3(x)
		hf <- function(x) hf3(x)
		gfd <- function(x) gfd3(x)
		}
	else if (link==4) {
		gf <- function(x) gf4(x)
		hf <- function(x) hf4(x)
		gfd <- function(x) gfd4(x)
		}
	else {cat("Link function is not defined")}
	ww <- function(mu, Sigma) {1/(Sigma*gfd(mu)^2)*(3*Sigma/pp(mu)+1/pp(mu)^3)}
    tf <- function(y, mu, Sigma) {gf(mu)+uu(y, mu)/(Sigma*ww(mu, Sigma)*gfd(mu))}
	norm <- function(x) {sqrt(sum(x^2))}

        y <- as.vector(y)
        X <- as.matrix(X)
        N <- length(y)
        p <- dim(X)[2] 
        q <- dim(Z)[2]
## initial values of beta, mu:
		if (length(beta) != p){
			mod1 <- glm(gf(y)~ X[, -1], family = gaussian)
			beta <- as.vector(mod1$coef)   
   	   	}	
   	   	mu <- as.vector(hf(X %*% beta))	
## initial values of alpha:
	    ddd <- dd(y, mu)
   	   	options(warn = -1)
	    mod2 <- try(glm(ddd ~ Z[,-1], family = Gamma(link = log)), TRUE)
   	   	if(inherits(mod2, "try-error"))
   	   	   	mod2 <- glm(ddd ~ Z[,-1], family = Gamma(link = log), start = c(mean(log(ddd)), 
				rep(0, q-1)), maxit = 5000)
   	   	options(warn = 0)
   	   	if (length(alpha) != q)	
   	   	   	alpha <- as.numeric(mod2$coef)
	    Sigma <- as.vector(exp(Z %*% alpha))
## Newton-Scoring algorithm for beta:
#        cat("initial beta-alpha:", round(c(beta, alpha), 4), "\n")
        iter <- 0
        dif <- 1
        beta.old <- beta
        while(dif > tol && iter < maxit) {
                iter <- iter + 1
#                cat("iter=", iter, "@", date(), "\n")
## 	Calculate Sigma and alpha:	
		         VF <- vf(y, mu, Sigma) 
                SF <- matrix(0, q, q)
                PhF <- rep(0, q)
                for(k in 1:N) { 
                    FV <- as.vector(Z[k, ])
                    SF <- SF + 0.5 * FV %*% t(FV)
                    PhF <- PhF + FV * VF[k]
                }
                invSF <- INV(SF)
                alpha <- as.vector(alpha + invSF %*% PhF)
                Sigma <- as.vector(exp(Z %*% alpha))
##       
                Y <- tf(y, mu, Sigma)
                eta <- gf(mu)
                w <- ww(mu, Sigma)
## Calculate Sensitivity matrix and scoring vector
                S1 <- matrix(0, p, p)
                Psi <- rep(0, p)
                for(k in 1:N) { 
                    vv <- as.vector(X[k, ])
                    S1 <- S1 + vv %*% as.matrix(w[k]) %*% t(vv)
                    Psi <- Psi + vv * as.numeric(w[k] * (Y[k] - eta[k]))
                }
                invS1 <- INV(S1)
                beta <- as.vector(beta.old + invS1 %*% Psi)
                mu <- as.vector(hf(X %*% beta))
                diff <- as.vector(beta - beta.old)
                dif <- norm(diff)/norm(beta)
#                cat("beta-alpha-dif:", round(c(beta, alpha, dif), 6), "\n")
                beta.old <- beta
        }
        sdfixefb <- sqrt(diag(invS1))
        fixefb <- as.matrix(cbind(beta, sdfixefb))
        VF <- vf(y, mu, Sigma) 
        SF <- matrix(0, q, q)
        PhF <- rep(0, q)
        for(k in 1:N) { 
             FV <- as.vector(Z[k, ])
             SF <- SF + 0.5 * FV %*% t(FV)
             PhF <- PhF + FV * VF[k]
        }
        invSF <- INV(SF)
        alpha <- as.vector(alpha + invSF %*% PhF)
        stddispar <- sqrt(diag(invSF))  
        fixefa <- as.matrix(cbind(alpha, stddispar))
        Sigma <- as.vector(exp(Z %*% alpha))       
##        out <- list(Dispersion = fixefa, Fixef = fixefb, Sigma=Sigma) 
##
		devi <- sum(dd(y,mu)/Sigma)
        ee <- as.vector((y - mu)/sqrt(mu * (1 - mu)))
        eee <- as.vector((y - mu)/sqrt(var.sim(mu, Sigma)))
        ss <- Sigma * ((3 * Sigma)/(mu * (1 - mu)) + 1/Vfun(mu))
        sss <- as.vector(gf(mu) + uu(y, mu)/sqrt(ss))
   	   	rr <- as.vector(uu(y, mu)/sqrt(Sigma * halfEd2(mu, Sigma)))
        pred <- as.vector(X %*% beta)
        loglike <- -sum(dd(y,mu)/Sigma - log(Sigma))/2
   	   	if (iter == maxit)
   	   	   	warning("step size truncated due to divergence")
#	Waldsigma <- t(sigma[2:q]) %*% INV(invGinf[(p + 2):(p + q), (p + 2):(p + q)]) %*% sigma[2:q]
   	return(list(fixef = fixefb, dispar = fixefa, Dispersion = Sigma, appstdPerr = ee, 
		stdPerr = eee, meanmu = mu, stdscor = rr, adjvar = sss, predict = pred, 
   	   	loglike = loglike, deviance = devi, iter = iter))
}
