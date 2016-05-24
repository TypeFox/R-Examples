.newt <- function(fn,start,...,eps.p = 1e-08,eps.v = NULL,
                  maxit = 50,verb = FALSE)
{
p.o <- start
itno <- 1
repeat {
	fj <- fn(p.o,...)
	v <- fj$fval
	t1 <- if(is.null(eps.v)) NULL else sum(abs(v))
	J <- as.matrix(fj$jacobian)
	if(qr(J)$rank < ncol(J)) {
		cat("Singular Jacobian.\n")
		rslt <- if(is.null(eps.v)) NA else if(t1 < eps.v) p.o
		else NA
		break
	}
	else {
		p.n <- p.o - solve(J) %*% v
		t2 <- max(abs(p.n - p.o))
		if(verb) {
			tmp <- format(round(c(p.o,p.n,v,t2,t1),6))
			np <- length(v)
			v1 <- paste(tmp[1:np],collapse = "  ")
			v2 <- paste(tmp[(np + 1):(2 * np)],collapse = "  ")
			v3 <- paste(tmp[(2 * np + 1):(3 * np)],collapse = "  ")
			v4 <- tmp[3 * np + 1]
			v5 <- tmp[3 * np + 2]
			cat("\nIteration  : ",itno,"\n",sep = "")
			cat("Old par    : ",v1,"\n",sep = "")
			cat("New par    : ",v2,"\n",sep = "")
			cat("Test ch.par: ",v4,"\n",sep = "")
			cat("Fn. vals.  : ",v3,"\n",sep = "")
			if(!is.null(t1))
			  cat("Test f.val: ",v5,"\n",sep = "")
		}
		if((!is.null(t1) && t1 < eps.v) | t2 < eps.p) {
			rslt <- p.n
			break
		}
		itno <- itno + 1
		if(itno > maxit) {
			cat("Newton's method failed to converge in\n")
			cat(maxit,"iterations.\n")
			rslt <- NA
			break
		}
		p.o <- p.n
	}
}
as.vector(rslt)
}
