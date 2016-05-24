metahet.base <- function(y, s2){
	if(length(y) != length(s2) | any(s2 < 0)) stop("error in the input data.")
	n <- length(y)

	w <- 1/s2
	mu.bar <- sum(w*y)/sum(w)

	out <- NULL

	out$weighted.mean <- mu.bar

	# the conventional methods
	Q <- sum(w*(y - mu.bar)^2)
	H <- sqrt(Q/(n - 1))
	I2 <- (Q - n + 1)/Q
	tau2.DL <- (Q - n + 1)/(sum(w) - sum(w^2)/sum(w))
	tau2.DL <- max(c(0, tau2.DL))
	out$Q <- Q
	out$H <- H
	out$I2 <- I2
	out$tau2.DL <- tau2.DL

	# absolute deviation based on weighted mean
	Qr <- sum(sqrt(w)*abs(y - mu.bar))
	Hr <- sqrt((3.14159*Qr^2)/(2*n*(n - 1)))
	Ir2 <- (Qr^2 - 2*n*(n - 1)/3.14159)/(Qr^2)
	tau2.r <- tau2.r.solver(w, Qr)
	out$Qr <- Qr
	out$Hr <- Hr
	out$Ir2 <- Ir2
	out$tau2.r <- tau2.r

	# absolute deviation based on weighted median

	expit <- function(x) {ifelse(x >= 0, 1/(1 + exp(-x/0.0001)), exp(x/0.0001)/(1 + exp(x/0.0001)))}
	psi <- function(x) {sum(w*(expit(x - y) - 0.5))}
	mu.med <- uniroot(psi, c(min(y) - 0.001, max(y) + 0.001))$root
	out$weighted.median <- mu.med
	Qm <- sum(sqrt(w)*abs(y - mu.med))
	Hm <- sqrt(3.14159/2)*Qm/n
	Im2 <- (Qm^2 - 2*n^2/3.14159)/Qm^2
	tau2.m <- tau2.m.solver(w, Qm)
	out$Qm <- Qm
	out$Hm <- Hm
	out$Im2 <- Im2
	out$tau2.m <- tau2.m

	return(out)
}