simglm1 <-
function(y, X, link=1, beta, maxit, tol) {
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
## initial values of beta, mu and Sigma:
	model <- glm(gf(y)~ X[, -1], family = gaussian)
   	if (length(beta) != p)
   	   	beta <- as.vector(model$coef)
	mu <- as.vector(hf(X %*% beta))
	Sigma <- sum(dd(y, mu))/(N - p)
## Newton-Scoring algorithm for beta:
	iter <- 0
	dif <- 1
	beta.old <- beta
	while(dif > tol && iter < maxit) {
		iter <- iter + 1
		Y <- tf(y, mu, Sigma)
		eta <- gf(mu)
		w <- ww(mu, Sigma)
		S1 <- matrix(0, p, p)
		Psi <- rep(0, p)
		for(k in 1:N) {
			vv <- as.vector(X[k,  ])
			S1 <- S1 + vv %*% as.matrix(w[k]) %*% t(vv)
			Psi <- Psi + vv * as.numeric(w[k] * (Y[k] - eta[k]))
		}
		beta <- as.vector(beta.old + INV(S1) %*% Psi)
		mu <- as.vector(hf(X %*% beta))
		Sigma <- sum(dd(y, mu))/(N - p)
		diff <- as.vector(beta - beta.old)
		dif <- norm(diff)/norm(beta)
		beta.old <- beta
	}
   	if (iter == maxit)
   	   	warning("step size truncated due to divergence")
	devi <- sum(dd(y,mu)/Sigma)
	sdbeta <- INV(S1)
	stdbeta <- sqrt(diag(sdbeta))
    fixefb <- cbind(beta, stdbeta)
    ee <- as.vector((y - mu)/sqrt(mu * (1 - mu)))
    eee <- as.vector((y - mu)/sqrt(var.sim(mu, Sigma)))
    ss <- Sigma * ((3 * Sigma)/(mu * (1 - mu)) + 1/Vfun(mu))
    sss <- as.vector(gf(mu) + uu(y, mu)/sqrt(ss))
    pred <- as.vector(X %*% beta)
	rr <- as.vector(uu(y, mu) / sqrt(Sigma * halfEd2(mu, rep(Sigma, N))))
    loglike <- -sum(dd(y,mu)/Sigma - log(Sigma))/2
   	return(list(fixef = fixefb, Dispersion = Sigma, appstdPerr = ee, stdPerr = eee, 
   	   	stdscor = rr, meanmu = mu, adjvar = sss, predict = pred, loglike = loglike, 
   	deviance = devi, iter = iter))
}
