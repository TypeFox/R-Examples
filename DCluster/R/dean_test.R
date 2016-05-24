# Copyright 2004 Virgilio Gómez Rubio and Roger Bivand

test.nb.pois <- function(x.nb, x.glm){
	if (!(inherits(x.nb, "negbin"))) stop("not a negbin object")
	if (!(inherits(x.glm, "glm"))) stop("not a glm object")
	zscore <- -log(x.nb$theta)*x.nb$theta/x.nb$SE.theta
	u <- pnorm(zscore)
	pz <- 2*min(u, 1-u)
	probs <- x.nb$theta/(x.nb$theta+x.nb$fitted.values)
	lrt <- 2*(sum(dnbinom(x.nb$y, x.nb$theta, probs, log=TRUE))-
		sum(dpois(x.nb$y, x.glm$fitted.values, log=TRUE)))
	names(lrt) <- "LR"
	pchi <- pchisq(lrt, df=1, lower.tail=FALSE)
	vec <- c(zscore, pz)
	names(vec) <- c("zscore", "p.mayor.modZ")
	res <- list(estimate=vec, statistic=lrt, p.value=pchi, parameter=1, 
		method="Likelihood ratio test for overdispersion",
		data.name=paste(deparse(substitute(x.nb)), ":", 
		deparse(substitute(x.glm))))
	class(res) <- "htest"
	res
}

DeanB <- function(x.glm, alternative="greater") {
	alternative <- match.arg(alternative, c("less", "greater", "two.sided"))
	if (!(inherits(x.glm, "glm"))) stop("not a glm object")
	y <- model.response(model.frame(x.glm))
	mu <- fitted(x.glm)
	Pb <- sum((y-mu)^2-y)/sqrt(2*sum(mu^2))
	names(Pb) <- "P_B"
	Pv <- NA
	if (is.finite(Pb)) {
			if (alternative == "two.sided") 
				Pv <- 2 * pnorm(abs(Pb), lower.tail=FALSE)
			else if (alternative == "greater")
				Pv <- pnorm(Pb, lower.tail=FALSE)
			else Pv <- pnorm(Pb)
	}
	res <- list(statistic=Pb, p.value=Pv, alternative=alternative, 
		method="Dean's P_B test for overdispersion",
		data.name=deparse(substitute(x.glm)))
	class(res) <- "htest"
	res
}

DeanB2 <- function(x.glm, alternative="greater"){
	
	y <- model.response(model.frame(x.glm))
	mu <- fitted(x.glm)
	h <- hatvalues(x.glm)
	Pb2 <- sum((y-mu)^2-y+h*mu)/sqrt(2*sum(mu^2))
	names(Pb2) <- "P'_B"
	Pv <- NA
	if (is.finite(Pb2)) {
			if (alternative == "two.sided") 
				Pv <- 2 * pnorm(abs(Pb2), lower.tail=FALSE)
			else if (alternative == "greater")
				Pv <- pnorm(Pb2, lower.tail=FALSE)
			else Pv <- pnorm(Pb2)
	}
	res <- list(statistic=Pb2, p.value=Pv, alternative=alternative, 
		method="Dean's P'_B test for overdispersion",
		data.name=deparse(substitute(x.glm)))
	class(res) <- "htest"
	res
}
