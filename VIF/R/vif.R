vif <- function(y, x, w0 = 0.05, dw = 0.05, subsize = 200, trace = TRUE, mode = c("dense", "sparse")){

    mf <- match.call();
    m <- match("mode", names(mf), 0L);
    if(m == 0)mode <- "sparse";

	y <- as.vector(y) - mean(y);
	n <- length(y);
	m <- subsize;
	if(m > n){
		print("m should be less than or equal to n"); 
		return(0);
	}
	x <- as.matrix(x);
	x <- t(t(x) - apply(x,2,mean));
	p <- dim(x)[2];
	sel <- rep(NA, 0);
	
	# initials
	res <- y - mean(y);
	sigma <- sd(res);
	w <- w0;
	flag <- 0;
	x.sel <- cbind(rep(1, n), x[, sel]);

	for(i in 1:p){
		x.can <- x[, i];
		if(mode == "dense")alpha <- w/(1 + i - flag) else alpha <- w/2/i;
		if(m == n)sam <- 1:n else sam <- sample(1:n, m);
		lm.can <- lm(x[sam, i] ~ x.sel[sam, ]-1);
#		if(!all(!is.na(lm.can$coef)))next;
		r2.x <- summary(lm.can)$r.squared;
		rho <- sqrt(1 - r2.x);
		t.abs <- abs(t(x.can) %*% res) / sqrt(t(x.can) %*% x.can) / rho / sigma;
		if(is.na(t.abs))next;
		pval <- 2 - 2 * pnorm(t.abs);
		if(trace == TRUE)print(paste(i, w, alpha, t.abs, pval, sep=" "));
		if(pval < alpha){
			x.tmp <- cbind(x.sel, x.can);
			lm.sel <- lm(y ~ x.tmp - 1);
			if(!all(!is.na(lm.sel$coef)))next;
			sel <- c(sel, i);
			x.sel <- x.tmp;
			res <- summary(lm.sel)$residuals;
			sigma <- sd(res) * sqrt(n - 1) / sqrt(n - 1 - length(sel));
			flag <- i;
			if(mode == "dense")w <- w + dw else w <- w + dw - alpha;
		} else if(mode == "dense")w <- w - alpha / (1 - alpha) else  w <- w - alpha;
		if(w <= 0)break;
	}
	
	return(list(select = sel, modelmatrix = x[, sel]));

}
