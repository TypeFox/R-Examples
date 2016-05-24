kInflatedLogConDiscr <- function(x, k = 0, prec1 = 1e-10, prec2 = 1e-15, itermax = 200, output = TRUE, theta0 = 0.5, p0 = NA){

	n <- length(x)
	out	<- sort(unique(x))
	z <- min(out):max(out)
    m <- length(z)	

	if((k %in% z) == FALSE){stop("k must be in the set {x_1, x1 + 1, ..., x_n - 1, x_n}")}
    if (identical(p0, NA)){p0 <- rep(1, diff(range(x)) + 1) / (diff(range(x)) + 1)}
    
    ## empirical on x1, x1 + 1, ..., xn - 1, xn
	emp <- rep(0, length(z))
	emp[out - min(x) + 1] <- as.vector(table(x) / n)

	theta.old <- theta0
	p.old <- p0
	w <- emp
    Eloglik.old <- prec1 + 1

	for(i in 1:itermax){

		theta.new <- emp[k - min(x) + 1] * theta.old / (theta.old + (1 - theta.old) * p.old[k - min(x) + 1]) 
		w[k - min(x) + 1] <- emp[k - min(x) + 1] - theta.new  

        # protect logConDiscrMLE() against instabilities due to zero values in w at the ends
		loc <- min(which(w > prec2)):max(which(w > prec2)) 
 
		mle <- logConDiscrMLE(z[loc], w = w[loc], output = FALSE)
		p.new <- rep(0, length(z))
		#p.new[mle$x - min(x) + 1] <- exp(mle$psiSupp)
        p.new[min(mle$x):max(mle$x) - min(x) + 1] <- exp(mle$psiSupp)
		Eloglik.new <- log(theta.new) * theta.new + log(1 - theta.new) * (1 - theta.new) + sum(w[p.new > 0] * log(p.new[p.new > 0]))
        Eloglik.diff <- Eloglik.new - Eloglik.old
        
      	if (identical(output, TRUE)){print(paste("iteration: ", format(i, digits = ceiling(log(itermax) / log(10))), " / E[comploglik]: ", format(Eloglik.new, digits = 4), " / diff loglik: ", format(Eloglik.diff, digits = 4), " / theta: ", format(theta.new, digits = 4), sep = ""))}
		if(abs(Eloglik.diff) < prec1){break}

		theta.old <- theta.new
		p.old <- p.new 
		Eloglik.old <- Eloglik.new
}

	q <- rep(0, m)
	q[k - min(x) + 1] <- 1
	f.fit <- theta.new * q + (1 - theta.new) * p.new
	res <- list("z" = z, "f" = f.fit, "E(L)" = Eloglik.new, "loglik" = sum(emp[f.fit > 0] * log(f.fit[f.fit > 0])), "theta" = theta.new, "logconc.pmf" = p.new, "logconc.z" = mle$x)
	return(res)
}


