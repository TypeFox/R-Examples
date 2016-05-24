max_threshold <- function(sample1, sample2, sig_mat){
	n1 <- dim(sample1)[1]; n2 <- dim(sample2)[1]; p <- dim(sample1)[2]
	a_f <- sqrt(2*log(log(p)))
	b_f <- 2*log(log(p)) + 0.5*log(log(log(p))) - 0.5*log(4*3.1416/(1 - 0.05)^2)
	eta <- 0.05
	diag.sig <- diag(sig_mat)
	diag.sig[diag.sig <= 10^(-10)] <- 10^(-10)
	T_orig <- (colMeans(sample1) - colMeans(sample2))^2/((1/n1 + 1/n2)*diag.sig)
	s_level <- T_orig[sign(T_orig) >= 0.01 & T_orig <= 2*(1 - eta)*log(p)]
	s_m <- matrix(s_level, length(s_level), p, byrow = F)
	T_m <- matrix((T_orig-1), length(s_level), p, byrow = T)
	T_m[sign(T_m + 1 - s_m) == -1] <- 0
	thr <- rowSums(T_m) 
	mean_thr <- 2*sqrt(s_level)*dnorm(sqrt(s_level))*p
	sd_thr <- sqrt(p*(2*((sqrt(s_level))^3 + sqrt(s_level))*dnorm(sqrt(s_level)) + 4 - 4*pnorm(sqrt(s_level))) - mean_thr^2/p)
	max_threshold <- max((thr-mean_thr)/sd_thr)*a_f - b_f
	return(max_threshold)            	
}
