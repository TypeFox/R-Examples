bayesnormtol.int <- function(x = NULL, norm.stats = list(x.bar=NA,s=NA,n=NA), alpha = 0.05, P = 0.99, side = 1, method = c("HE", "HE2", "WBE", "ELL", "KM", 
            "EXACT", "OCT"), m = 50, hyper.par = list(mu.0=NULL,sig2.0=NULL,m.0=NULL,n.0=NULL)){
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
	method <- match.arg(method)
	if (is.null(x)){
		x.bar <- norm.stats$x.bar
		s <- norm.stats$s
		n <- norm.stats$n
	} else{
		x.bar <- mean(x)
		s <- sd(x)
		n <- length(x)	
	}
		if (sum(is.na(hyper.par))%in%c(1:3)) {
	        stop(paste("All hyperparameters must be specified!", 
	            "\n"))
		} else if(all(is.na(hyper.par))) { 
			    K <- invisible(K.factor(n = n, alpha = alpha, P = P, side = side, method = method, m = m))
				lower <- x.bar - s*K
				upper <- x.bar + s*K
		} else { 
				mu.0 <- hyper.par$mu.0; sig2.0 <- hyper.par$sig2.0; m.0 <- hyper.par$m.0; n.0 <- hyper.par$n.0
			    K <- invisible(K.factor(n = n.0 + n, f = m.0 + n - 1, alpha = alpha, P = P, side = side, method = method, m = m))
				x.bar.bar <- (n.0*mu.0+n*x.bar)/(n.0+n)
				q2 <- (m.0*sig2.0+(n-1)*s^2+(n.0*n*(x.bar-mu.0)^2)/(n.0+n))/(m.0+n-1)
				lower <- x.bar.bar - sqrt(q2)*K
				upper <- x.bar.bar + sqrt(q2)*K
		}
    if (side == 1) {
        temp <- data.frame(cbind(alpha, P, lower, upper))
        colnames(temp) <- c("alpha", "P", "1-sided.lower", 
            "1-sided.upper")
    }
    else {
        temp <- data.frame(cbind(alpha, P, lower, upper))
        colnames(temp) <- c("alpha", "P", "2-sided.lower", 
            "2-sided.upper")
    }
    temp
}
