# last modified 21 Oct 08 by J. Fox

"polyserial" <-
	function(x, y, ML=FALSE, control=list(), std.err=FALSE, maxcor=.9999, bins=4){
	f <- function(pars){
		rho <- pars[1]
		if (abs(rho) > maxcor) rho <- sign(rho)*maxcor
		cts <- if (length(pars) == 1) c(-Inf, cuts, Inf)
			else c(-Inf, pars[-1], Inf)
		tau <- (matrix(cts, n, s+1, byrow=TRUE) - matrix(rho*z, n, s+1))/
			sqrt(1 - rho^2)
		- sum(log(dnorm(z)*(pnorm(tau[cbind(indices, y+1)]) - pnorm(tau[cbind(indices, y)]))))
	}
	if (!is.numeric(x)) stop("x must be numeric")
	valid <- complete.cases(x, y)
	x <- x[valid]
	y <- y[valid]
	z <- scale(x)
	tab <- table(y)
	n <- sum(tab)
	s <- length(tab)
	if (s < 2) {
		warning("y has fewer than 2 levels")
		return(NA)
	}
	if (sum(0 != tab) < 2){
		warning("y has fewer than 2 levels with data")
		return(NA)
	}
	indices <- 1:n
	cuts <- qnorm(cumsum(tab)/n)[-s]
	y <- as.numeric(as.factor(y))
	rho <- sqrt((n - 1)/n)*sd(y)*cor(x, y)/sum(dnorm(cuts))
	if (ML) {
		result <- optim(c(rho, cuts), f, control=control, hessian=std.err)
		if (result$par[1] > 1){
			result$par[1] <- 1
			warning("inadmissible correlation set to 1")
		}
		else if (result$par[1] < -1){
			result$par[1] <- -1
			warning("inadmissible correlation set to -1")
		}
		if (std.err){
			chisq <- chisq(y, z, result$par[1], result$par[-1], bins=bins)
			df <- s*bins - s  - 1
			result <- list(type="polyserial",
				rho=result$par[1],
				cuts=result$par[-1],
				var=solve(result$hessian),
				n=n,
				chisq=chisq,
				df=df,
				ML=TRUE)
			class(result) <- "polycor"
			return(result)
		}
		else return(as.vector(result$par[1]))  
	}
	else if (std.err){
		result <- optim(rho, f, control=control, hessian=TRUE, method="BFGS")
		if (result$par > 1){
			result$par <- 1
			warning("inadmissible correlation set to 1")
		}
		else if (result$par < -1){
			result$par <- -1
			warning("inadmissible correlation set to -1")
		}
		chisq <- chisq(y, z, rho, cuts, bins=bins)
		df <- s*bins - s  - 1
		result <- list(type="polyserial",
			rho=result$par,
			var=1/result$hessian,
			n=n,
			chisq=chisq,
			df=df,
			ML=FALSE)
		class(result) <- "polycor"
		return(result)
	}
	else {
		if (rho > 1){
			rho <- 1
			warning("inadmissible correlation set to 1")
		}
		else if (rho < -1){
			rho <- -1
			warning("inadmissible correlation set to -1")
		}
		rho
	}
}
