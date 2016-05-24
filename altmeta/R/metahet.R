metahet <- function(y, s2, n.resam = 1000){
	if(length(y) != length(s2) | any(s2 < 0)) stop("error in the input data.")
	n <- length(y)

	out0 <- metahet.base(y = y, s2 = s2)

	out.resam.param <- lapply(1:n.resam, function(o){
		y.temp <- rnorm(n, mean = out0$weighted.mean, sd = sqrt(s2))
		out.temp <- metahet.base(y.temp, s2)
		return(out.temp)
	})
	out.resam.param <- unlist(out.resam.param)

	out.resam.replace <- lapply(1:n.resam, function(o){
		temp <- sample(n, replace = TRUE)
		y.temp <- y[temp]
		s2.temp <- s2[temp]
		out.temp <- metahet.base(y.temp, s2.temp)
		return(out.temp)
	})
	out.resam.replace <- unlist(out.resam.replace)

	# the conventional methods
	Q <- out0$Q
	Q.resam.param <- out.resam.param[names(out.resam.param) == "Q"]
	names(Q.resam.param) <- NULL
	p.Q <- (sum(Q.resam.param > Q) + 1)/(n.resam + 1)
	p.Q.theo <- 1 - pchisq(Q, n - 1)

	Q.resam.replace <- out.resam.replace[names(out.resam.replace) == "Q"]
	names(Q.resam.replace) <- NULL
	ci.Q <- quantile(Q.resam.replace, probs = c(0.025, 0.975))

	tau2.DL <- out0$tau2.DL
	tau2.resam.replace <- out.resam.replace[names(out.resam.replace) == "tau2.DL"]
	names(tau2.resam.replace) <- NULL
	ci.tau2.DL <- quantile(tau2.resam.replace, probs = c(0.025, 0.975))

	H <- out0$H
	H.resam.replace <- out.resam.replace[names(out.resam.replace) == "H"]
	names(H.resam.replace) <- NULL
	ci.H <- quantile(H.resam.replace, probs = c(0.025, 0.975))

	I2 <- out0$I2
	I2.resam.replace <- out.resam.replace[names(out.resam.replace) == "I2"]
	names(I2.resam.replace) <- NULL
	ci.I2 <- quantile(I2.resam.replace, probs = c(0.025, 0.975))

	# absolute deviation based on weighted mean
	Qr <- out0$Qr
	Qr.resam.param <- out.resam.param[names(out.resam.param) == "Qr"]
	names(Qr.resam.param) <- NULL
	p.Qr <- (sum(Qr.resam.param > Qr) + 1)/(n.resam + 1)

	Qr.resam.replace <- out.resam.replace[names(out.resam.replace) == "Qr"]
	names(Qr.resam.replace) <- NULL
	ci.Qr <- quantile(Qr.resam.replace, probs = c(0.025, 0.975))

	tau2.r <- out0$tau2.r
	tau2.r.resam.replace <- out.resam.replace[names(out.resam.replace) == "tau2.r"]
	names(tau2.r.resam.replace) <- NULL
	ci.tau2.r <- quantile(tau2.r.resam.replace, probs = c(0.025, 0.975))

	Hr <- out0$Hr
	Hr.resam.replace <- out.resam.replace[names(out.resam.replace) == "Hr"]
	names(Hr.resam.replace) <- NULL
	ci.Hr <- quantile(Hr.resam.replace, probs = c(0.025, 0.975))

	Ir2 <- out0$Ir2
	Ir2.resam.replace <- out.resam.replace[names(out.resam.replace) == "Ir2"]
	names(Ir2.resam.replace) <- NULL
	ci.Ir2 <- quantile(Ir2.resam.replace, probs = c(0.025, 0.975))

	# absolute deviation based on weighted median
	Qm <- out0$Qm
	Qm.resam.param <- out.resam.param[names(out.resam.param) == "Qm"]
	names(Qm.resam.param) <- NULL
	p.Qm <- (sum(Qm.resam.param > Qm) + 1)/(n.resam + 1)

	Qm.resam.replace <- out.resam.replace[names(out.resam.replace) == "Qm"]
	names(Qm.resam.replace) <- NULL
	ci.Qm <- quantile(Qm.resam.replace, probs = c(0.025, 0.975))

	tau2.m <- out0$tau2.m
	tau2.m.resam.replace <- out.resam.replace[names(out.resam.replace) == "tau2.m"]
	names(tau2.m.resam.replace) <- NULL
	ci.tau2.m <- quantile(tau2.m.resam.replace, probs = c(0.025, 0.975))

	Hm <- out0$Hm
	Hm.resam.replace <- out.resam.replace[names(out.resam.replace) == "Hm"]
	names(Hm.resam.replace) <- NULL
	ci.Hm <- quantile(Hm.resam.replace, probs = c(0.025, 0.975))

	Im2 <- out0$Im2
	Im2.resam.replace <- out.resam.replace[names(out.resam.replace) == "Im2"]
	names(Im2.resam.replace) <- NULL
	ci.Im2 <- quantile(Im2.resam.replace, probs = c(0.025, 0.975))

	# output
	out <- NULL

	out$p.Q <- p.Q
	out$p.Q.theo <- p.Q.theo
	out$p.Qr <- p.Qr
	out$p.Qm <- p.Qm

	out$Q <- Q
	out$ci.Q <- ci.Q
	out$tau2.DL <- tau2.DL
	out$ci.tau2.DL <- ci.tau2.DL
	out$H <- H
	out$ci.H <- ci.H
	out$I2 <- I2
	out$ci.I2 <- ci.I2

	out$Qr <- Qr
	out$ci.Qr <- ci.Qr
	out$tau2.r <- tau2.r
	out$ci.tau2.r <- ci.tau2.r
	out$Hr <- Hr
	out$ci.Hr <- ci.Hr
	out$Ir2 <- Ir2
	out$ci.Ir2 <- ci.Ir2

	out$Qm <- Qm
	out$ci.Qm <- ci.Qm
	out$tau2.m <- tau2.m
	out$ci.tau2.m <- ci.tau2.m
	out$Hm <- Hm
	out$ci.Hm <- ci.Hm
	out$Im2 <- Im2
	out$ci.Im2 <- ci.Im2

	return(out)
}