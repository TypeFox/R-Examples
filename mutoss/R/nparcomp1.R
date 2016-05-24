# Simultaneous confidence intervals for relative contrast effects
# 
# Author: FrankKonietschke
###############################################################################


nparcomp <- function (formula, data, type = c("UserDefined", "Tukey", "Dunnett", "Sequen",
				"Williams", "Changepoint", "AVE", "McDermott", "Marcus","UmbrellaWilliams"),
		control = NULL, conflevel = 0.95, alternative = c("two.sided",
				"less", "greater"), rounds = 3, correlation = FALSE,
		asy.method = c("logit", "probit", "normal", "mult.t"), plot.simci = FALSE,
		info = TRUE, contrastMatrix = NULL)
{
	corr.mat <- function(m, nc) {
		rho <- matrix(c(0), ncol = nc, nrow = nc)
		for (i in 1:nc) {
			for (j in 1:nc) {
				rho[i, j] <- m[i, j]/sqrt(m[i, i] * m[j, j])
			}
		}
		return(rho)
	}
	ssq <- function(x) {
		sum(x * x)
	}
	logit <- function(p) {
		log(p/(1 - p))
	}
	probit <- function(p) {
		qnorm(p)
	}
	expit <- function(G) {
		exp(G)/(1 + exp(G))
	}
	index <- function(char, test) {
		nc <- length(char)
		for (i in 1:nc) {
			if (char[i] == test) {
				return(i)
			}
		}
	}
	z.quantile <- function(conflevel = conflevel, corr, a, df = df.sw,
			dbs) {
		if (dbs == "n") {
			if (a == "two.sided") {
				z <- qmvnorm(conflevel, corr = corr, tail = "both")$quantile
			}
			if (a == "less" || a == "greater") {
				z <- qmvnorm(conflevel, corr = corr, tail = "lower")$quantile
			}
		}
		if (dbs == "t") {
			if (a == "two.sided") {
				z <- qmvt(conflevel, df = df.sw, interval = c(-10,
								10), corr = corr, tail = "both")$quantile
			}
			if (a == "less" || a == "greater") {
				z <- qmvt(conflevel, df = df.sw, interval = c(-10,
								10), corr = corr, tail = "lower")$quantile
			}
		}
		return(z)
	}
	if (conflevel >= 1 || conflevel <= 0) {
		stop("The confidence level must be between 0 and 1!")
		if (is.null(alternative)) {
			stop("Please declare the alternative! (two.sided, less, greater)")
		}
	}
	type <- match.arg(type)
	alternative <- match.arg(alternative)
	asy.method <- match.arg(asy.method)
	if (length(formula) != 3) {
		stop("You can only analyse one-way layouts!")
	}
	dat <- model.frame(formula, data)
	if (ncol(dat) != 2) {
		stop("Specify one response and only one class variable in the formula")
	}
	if (is.numeric(dat[, 1]) == FALSE) {
		stop("Response variable must be numeric")
	}
	response <- dat[, 1]
	factorx <- as.factor(dat[, 2])
	fl <- levels(factorx)
	a <- nlevels(factorx)
	if (a <= 2) {
		stop("You want to perform a two-sample test. Please use the function npar.t.test")
	}
	samples <- split(response, factorx)
	n <- sapply(samples, length)
	if (any(n <= 1)) {
		warn <- paste("The factor level", fl[n <= 1], "has got only one observation!")
		stop(warn)
	}
	ntotal <- sum(n)
	a <- length(n)
	tmp <- expand.grid(1:a, 1:a)
	ind <- tmp[[1]] > tmp[[2]]
	vi <- tmp[[2]][ind]
	vj <- tmp[[1]][ind]
	nc <- length(vi)
	gn <- n[vi] + n[vj]
	intRanks <- lapply(samples, rank)
	pairRanks <- lapply(1:nc, function(arg) {
				rank(c(samples[[vi[arg]]], samples[[vj[arg]]]))
			})
	pd <- sapply(1:nc, function(arg) {
				i <- vi[arg]
				j <- vj[arg]
				(sum(pairRanks[[arg]][(n[i] + 1):gn[arg]])/n[j] - (n[j] +
								1)/2)/n[i]
			})
	dij <- dji <- list(0)
	sqij <- sapply(1:nc, function(arg) {
				i <- vi[arg]
				j <- vj[arg]
				pr <- pairRanks[[arg]][(n[i] + 1):gn[arg]]
				dij[[arg]] <<- pr - sum(pr)/n[j] - intRanks[[j]] + (n[j] +
							1)/2
				ssq(dij[[arg]])/(n[i] * n[i] * (n[j] - 1))
			})
	sqji <- sapply(1:nc, function(arg) {
				i <- vi[arg]
				j <- vj[arg]
				pr <- pairRanks[[arg]][1:n[i]]
				dji[[arg]] <<- pr - sum(pr)/n[i] - intRanks[[i]] + (n[i] +
							1)/2
				ssq(dji[[arg]])/(n[j] * n[j] * (n[i] - 1))
			})
	vd.bf <- ntotal * (sqij/n[vj] + sqji/n[vi])
	singular.bf <- (vd.bf == 0)
	vd.bf[singular.bf] <- 1e-05
	df.sw <- (n[vi] * sqij + n[vj] * sqji)^2/((n[vi] * sqij)^2/(n[vj] -
					1) + (n[vj] * sqji)^2/(n[vi] - 1))
	lambda <- sqrt(n[vi]/(gn + 1))
	cov.bf1 <- diag(nc)
	rho.bf <- diag(nc)
	for (x in 1:(nc - 1)) {
		for (y in (x + 1):nc) {
			i <- vi[x]
			j <- vj[x]
			v <- vi[y]
			w <- vj[y]
			p <- c(i == v, j == w, i == w, j == v)
			if (sum(p) == 1) {
				cl <- list(function() (t(dji[[x]]) %*% dji[[y]])/(n[j] *
										n[w] * n[i] * (n[i] - 1)), function() (t(dij[[x]]) %*%
										dij[[y]])/(n[i] * n[v] * n[j] * (n[j] - 1)),
						function() -(t(dji[[x]]) %*% dij[[y]])/(n[v] *
										n[j] * n[i] * (n[i] - 1)), function() -(t(dij[[x]]) %*%
										dji[[y]])/(n[i] * n[w] * n[j] * (n[j] - 1)))
				case <- (1:4)[p]
				rho.bf[x, y] <- rho.bf[y, x] <- sqrt(ntotal *
								ntotal)/sqrt(vd.bf[x] * vd.bf[y]) * cl[[case]]()
				cov.bf1[x, y] <- cov.bf1[y, x] <- sqrt(vd.bf[x] *
								vd.bf[y])
			}
		}
	}
	V <- (cov.bf1 + diag(vd.bf - 1)) * rho.bf
	cov.bf1 <- cbind(V,-1*V)
	cov.bf2<-cbind(-1*V,V)
	cov.bf <- rbind(cov.bf1,cov.bf2)
	switch(type,
			UserDefined = {
				if (is.null(contrastMatrix)){stop("Give a contrast matrix with the contrast.matrix = <> option or choose a contrast!")}
				
				if (is.null(control)) {
					nc <- nrow(contrastMatrix)
					weights.help <- weightMatrix(n, contrast.matrix = contrastMatrix)
					weight<-weights.help$weight.matrix
					weight.help <- weights.help$weight.help
					cmpid <- paste("C", 1:nc)
					type.of.contrast <- "User Defined"
				}
				else {
					stop("Please declare the control group via your contrast matrix!")
				}
				type.of.contrast <- "User Defined"
			},
			
			Tukey = {
				if (is.null(control)) {
					if (alternative!="two.sided"){stop("The Tukey contrast can only be tested two-sided!")}
					nc <- a * (a - 1)/2
					cmpid <- sapply(1:nc, function(arg) {
								i <- vi[arg]
								j <- vj[arg]
								paste("p", "(", fl[i], ",", fl[j], ")", sep = "")
							})
					
					weights.help <- weightMatrix(n, "Tukey")
					weight<-weights.help$weight.matrix
					weight.help <- weights.help$weight.help
				}
				else {
					stop("The Tukey contrast hasn't got a control group!")
				}
				type.of.contrast <- "Tukey"
			}, Dunnett = {
				nc <- a - 1
				if (is.null(control)) {
					cont <- 1
				}
				else {
					if (!any(fl == control)) {
						stop("The dataset doesn't contain this control group!")
					}
					cont <- which(fl == control)
				}
				vj <- which((1:a) != cont)
				vi <- rep(cont, a - 1)
				weights.help <- weightMatrix(n, "Dunnett",cont)
				weight<-weights.help$weight.matrix
				weight.help <- weights.help$weight.help
				cmpid <- sapply(1:nc, function(arg) {
							i <- vi[arg]
							j <- vj[arg]
							paste("p", "(", fl[i], ",", fl[j], ")", sep = "")
						})
				
				type.of.contrast <- "Dunnett"
			}, Sequen = {
				if (is.null(control)) {
					nc <- a - 1
					vi <- 1:(a - 1)
					vj <- 2:a
					weights.help <- weightMatrix(n, "Sequen")
					weight<-weights.help$weight.matrix
					weight.help <- weights.help$weight.help
					cmpid <- sapply(1:nc, function(arg) {
								i <- vj[arg]
								j <- vi[arg]
								paste("p", "(", fl[j], ",", fl[i], ")", sep = "")
							})
					
					type.of.contrast <- "Sequen"
				}
				else {
					stop("The Sequen-Contrast hasn't got a control group!")
				}
			},
			
			
			Williams = {
				if (is.null(control)) {
					
					nc <- a - 1
					weights.help <- weightMatrix(n, "Williams")
					weight<-weights.help$weight.matrix
					weight.help <- weights.help$weight.help
					cmpid <- paste("C", 1:(a - 1))
					type.of.contrast <- "Williams"
				}
				else {
					stop("The Williams contrast hasn't got a control group!")
				}
			}, Changepoint = {
				if (is.null(control)) {
					nc <- a - 1
					weights.help <- weightMatrix(n, "Changepoint")
					weight<-weights.help$weight.matrix
					weight.help <- weights.help$weight.help
					cmpid <- paste("C", 1:(a - 1))
					type.of.contrast <- "Changepoint"
				}
				else {
					stop("The Changepoint-Contrast hasn't got a control group!")
				}
			}, AVE = {
				if (is.null(control)) {
					nc <- a
					weights.help <- weightMatrix(n, "AVE")
					weight<-weights.help$weight.matrix
					weight.help <- weights.help$weight.help
					cmpid <- paste("C", 1:a)
					type.of.contrast <- "Average"
				}
				else {
					stop("The Average-Contrast hasn't got a control group!")
				}
			}, McDermott = {
				if (is.null(control)) {
					nc <- a - 1
					weights.help <- weightMatrix(n, "McDermott")
					weight<-weights.help$weight.matrix
					weight.help <- weights.help$weight.help
					cmpid <- paste("C", 1:(a - 1))
					
					type.of.contrast <- "McDermott"
				}
				else {
					stop("The McDermott-Contrast hasn't got a control group!")
				}
			}, UmbrellaWilliams = {
				if (is.null(control)) {
					nc <- a * (a - 1)/2
					weights.help <- weightMatrix(n, "UmbrellaWilliams")
					weight<-weights.help$weight.matrix
					weight.help <- weights.help$weight.help
					cmpid <- paste("C", 1:(a * (a - 1)/2))
					
					type.of.contrast <- "Umbrella Williams"
				}
				else {
					stop("The Umbrella Williams-Contrast hasn't got a control group!")
				}
			},
			
			Marcus = {
				if (is.null(control)) {
					nc <- a * (a - 1)/2
					weights.help <- weightMatrix(n, "Marcus")
					weight<-weights.help$weight.matrix
					weight.help <- weights.help$weight.help
					cmpid <- paste("C", 1:(a * (a - 1)/2))
					
					type.of.contrast <- "Marcus"
				}
				else {
					stop("The Marcus-Contrast hasn't got a control group!")
				}
			}
	
	
	
	)
	
	pd1 <- (pd == 1)
	pd0 <- (pd == 0)
	pd[pd1] <- 0.999
	pd[pd0] <- 0.001
	pd.help1 <- c(pd, 1-pd)
	pd <- c(weight%*%pd.help1)
	cov.bf <- weight %*% cov.bf %*% t(weight)
	for (i in 1:nc) {
		if (cov.bf[i, i] == 0) {
			cov.bf[i, i] <- 0.001
		}
	}
	
	vd.bf <- c(diag(cov.bf))
	vd.bf <- c(vd.bf)
	rho.bf <- corr.mat(cov.bf, nc)
	t.bf <- sqrt(ntotal) * (pd - 1/2)/sqrt(vd.bf)
	rownames(weight) <- paste("C", 1:nc)
	ncomp <- a * (a - 1)/2
	tmp <- expand.grid(1:a, 1:a)
	ind <- tmp[[1]] > tmp[[2]]
	v2 <- tmp[[2]][ind]
	v1 <- tmp[[1]][ind]
	namen1 <- sapply(1:ncomp, function(arg) {
				i <- v2[arg]
				j <- v1[arg]
				paste("p", "(", fl[i], ",", fl[j], ")", sep = "")
			})
	namen2 <- sapply(1:ncomp, function(arg) {
				i <- v2[arg]
				j <- v1[arg]
				paste("p", "(", fl[j], ",", fl[i], ")", sep = "")
			})
	colnames(weight) <- c(namen1,namen2)
	df.sw[is.nan(df.sw)] <- 1000
	df.sw <- weight %*% c(df.sw,df.sw)
	df.sw <- max(4, min(df.sw))
	pd <- c(pd)
	logit.pd <- logit(pd)
	logit.dev <- diag(1/(pd * (1 - pd)))
	logit.cov <- logit.dev %*% cov.bf %*% t(logit.dev)
	vd.logit <- c(diag(logit.cov))
	t.logit <- (logit.pd) * sqrt(ntotal/vd.logit)
	probit.pd <- qnorm(pd)
	probit.dev <- diag(sqrt(2 * pi)/(exp(-0.5 * qnorm(pd) * qnorm(pd))))
	probit.cov <- probit.dev %*% cov.bf %*% t(probit.dev)
	vd.probit <- c(diag(probit.cov))
	t.probit <- (probit.pd) * sqrt(ntotal/vd.probit)
	p.bfn = p.bft = p.bflogit = p.bfprobit = c()
	p.n = p.t = p.logit = p.probit = c()
	
	if (alternative == "two.sided") {
		z.bft <- z.quantile(conflevel = conflevel, corr = rho.bf,
				"two.sided", df = df.sw, dbs = "t")
		z.bfn <- z.quantile(conflevel = conflevel, corr = rho.bf,
				"two.sided", df = 0, dbs = "n")
		lower.bft <- pd - sqrt(vd.bf/ntotal) * z.bft
		upper.bft <- pd + sqrt(vd.bf/ntotal) * z.bft
		lower.bfn <- pd - sqrt(vd.bf/ntotal) * z.bfn
		upper.bfn <- pd + sqrt(vd.bf/ntotal) * z.bfn
		lower.logit <- expit(logit.pd - sqrt(vd.logit/ntotal) *
						z.bfn)
		upper.logit <- expit(logit.pd + sqrt(vd.logit/ntotal) *
						z.bfn)
		lower.probit <- pnorm(probit.pd - sqrt(vd.probit/ntotal) *
						z.bfn)
		upper.probit <- pnorm(probit.pd + sqrt(vd.probit/ntotal) *
						z.bfn)
		#---------------------------------------------------------#
		#--------------Berechnung von p-Werten -------------------#
		for (i in 1:nc) {
			p.bft[i] <- 1-pmvt(lower=-abs(t.bf[i]), abs(t.bf[i]), df=df.sw, corr=rho.bf,delta=rep(0,nc))
			p.bfn[i] <- 1-pmvnorm(lower=-abs(t.bf[i]), abs(t.bf[i]),corr=rho.bf,mean=rep(0,nc))
			p.bflogit[i] <- 1-pmvnorm(lower=-abs(t.logit[i]), abs(t.logit[i]),corr=rho.bf ,mean=rep(0,nc))
			p.bfprobit[i] <- 1-pmvnorm(lower=-abs(t.probit[i]), abs(t.probit[i]),corr=rho.bf,mean=rep(0,nc))
			p.tt <- pt(t.bf[i], df.sw)
			p.t[i] <- min(2*p.tt,2-2*p.tt)
			p.nn <- pnorm(t.bf[i])
			p.n[i] <- min(2*p.nn, 2-2*p.nn)
			p.ll <- pnorm(t.logit[i])
			p.logit[i] <- min(2*p.ll, 2-2*p.ll)
			p.pp <- pnorm(t.probit[i])
			p.probit[i] <- min(2*p.pp, 2-2*p.pp)
		}
		
		text.output.p <- "H_0: p(i,j)=1/2"
		text.output.KI <- paste(100 * conflevel, "%", "2-sided",
				"Simultaneous-Confidence-Intervals for Relative Effects")
		upper <- "]"
		lower <- "["
	}
	if (alternative == "greater") {
		z.bft <- z.quantile(conflevel = conflevel, corr = rho.bf,
				"less", df = df.sw, dbs = "t")
		z.bfn <- qmvnorm(conflevel, corr = rho.bf, tail = "lower")$quantile
		lower.bft <- pd - sqrt(vd.bf/ntotal) * z.bft
		lower.bfn <- pd - sqrt(vd.bf/ntotal) * z.bfn
		lower.logit <- expit(logit.pd - sqrt(vd.logit/ntotal) *
						z.bfn)
		lower.probit <- pnorm(probit.pd - sqrt(vd.probit/ntotal) *
						z.bfn)
		upper.bft = upper.probit = upper.logit = upper.bfn = 1
		for (i in 1:nc) {
			p.bfn[i] <- 1 - pmvnorm(lower = -Inf, upper = t.bf[i],
					mean = rep(0, nc), corr = rho.bf)
			p.bft[i] <- 1 - pmvt(lower = -Inf, upper = t.bf[i],
					delta = rep(0, nc), df = df.sw, corr = rho.bf)
			p.bflogit[i] <- 1 - pmvnorm(lower = -Inf, upper = t.logit[i],
					mean = rep(0, nc), corr = rho.bf)
			p.bfprobit[i] <- 1 - pmvnorm(lower = -Inf, upper = t.probit[i],
					mean = rep(0, nc), corr = rho.bf)
			
			p.t[i] <-1-pt(t.bf[i], df.sw)
			p.n[i] <- 1-pnorm(t.bf[i])
			p.logit[i] <- 1- pnorm(t.logit[i])
			p.probit[i] <-1- pnorm(t.probit[i])
			
		}
		text.output.p <- "H_0: p(i,j)<=1/2"
		text.output.KI <- paste(100 * conflevel, "%", "1-sided",
				"Simultaneous-Confidence-Intervals for Relative Effects")
		upper <- "]"
		lower <- "("
	}
	if (alternative == "less") {
		z.bft <- z.quantile(conflevel = conflevel, corr = rho.bf,
				"less", df = df.sw, dbs = "t")
		z.bfn <- z.quantile(conflevel = conflevel, corr = rho.bf,
				"less", df = 0, dbs = "n")
		upper.bft <- pd + sqrt(vd.bf/ntotal) * z.bft
		upper.bfn <- pd + sqrt(vd.bf/ntotal) * z.bfn
		upper.logit <- expit(logit.pd + sqrt(vd.logit/ntotal) *
						z.bfn)
		upper.probit <- pnorm(probit.pd + sqrt(vd.probit/ntotal) *
						z.bfn)
		lower.bft = lower.probit = lower.logit = lower.bfn = 0
		for (i in 1:nc) {
			p.bfn[i] <- 1 - pmvnorm(lower = t.bf[i], upper = Inf,
					mean = rep(0, nc), corr = rho.bf)
			p.bft[i] <- 1 - pmvt(lower = t.bf[i], upper = Inf,
					delta = rep(0, nc), df = df.sw, corr = rho.bf)
			p.bflogit[i] <- 1 - pmvnorm(lower = t.logit[i], upper = Inf,
					mean = rep(0, nc), corr = rho.bf)
			p.bfprobit[i] <- 1 - pmvnorm(lower = t.probit[i],
					upper = Inf, mean = rep(0, nc), corr = rho.bf)
			p.t[i] <-  pt(t.bf[i], df.sw)
			p.n[i] <- pnorm(t.bf[i])
			p.logit[i] <- pnorm(t.logit[i])
			p.probit[i] <- pnorm(t.probit[i])
		}
		text.output.p <- " H_0: p(i,j)>=1/2"
		text.output.KI <- paste(100 * conflevel, "%", "1-sided",
				"Simultaneous-Confidence-Intervals for Relative Effects")
		upper <- ")"
		lower <- "["
	}
	bfn.lower <-  round(lower.bfn, rounds)
	bfn.upper <- round(upper.bfn, rounds)
	bft.lower <-  round(lower.bft, rounds)
	bft.upper <- round(upper.bft,rounds)
	logit.lower <- round(lower.logit, rounds)
	logit.upper <- round(upper.logit,rounds)
	probit.lower <- round(lower.probit, rounds)
	probit.upper <- round(upper.probit,rounds)
	p.bflogit <- round(p.bflogit, rounds)
	p.bfprobit <- round(p.bfprobit, rounds)
	p.bft <- round(p.bft, rounds)
	p.bfn <- round(p.bfn, rounds)
	p.logit <- round(p.logit, rounds)
	p.probit <- round(p.probit, rounds)
	p.t <- round(p.t, rounds)
	p.n <- round(p.n, rounds)
	pd <- round(pd, rounds)
	if (correlation == TRUE) {
		Correlation <- list(Correlation.matrix.N = rho.bf, Covariance.matrix.N = cov.bf,
				Warning = paste("Attention! The covariance matrix is multiplied with N",
						"=", ntotal))
	}
	else {
		Correlation <- NA
	}
	data.info <- data.frame(row.names = 1:a, Sample = fl, Size = n)
	switch(asy.method, logit = {
				x.werte = cbind(lower.logit, pd, upper.logit)
				result <- list(weight.matrix = weight, Data.Info = data.info,
						Analysis.of.relative.effects = data.frame(row.names = c(1:nc),
								comparison = cmpid, rel.effect = pd, lower= logit.lower, upper=logit.upper,
								t.value = t.logit, p.adj = p.bflogit, p.raw= p.logit), Mult.Distribution = data.frame(Quantile = z.bfn,
								p.Value.global = min(p.bflogit)), Correlation = Correlation)
				Asymptotic.Method <- "Multivariate Delta-Method (Logit)"
			}, probit = {
				x.werte = cbind(lower.probit, pd, upper.probit)
				result <- list(weight.matrix = weight, Data.Info = data.info,
						Analysis.of.relative.effects = data.frame(row.names = c(1:nc),
								comparison = cmpid, rel.effect = pd, lower = probit.lower, upper = probit.upper,
								t.value = t.probit, p.adj = p.bfprobit, p.raw = p.probit), Mult.Distribution = data.frame(Quantile = z.bfn,
								p.Value.global = min(p.bfprobit)), Correlation = Correlation)
				Asymptotic.Method <- "Multivariate Delta-Method (Probit)"
			}, normal = {
				x.werte = cbind(lower.bfn, pd, upper.bfn)
				result <- list(weight.matrix = weight, Data.Info = data.info,
						Analysis.of.relative.effects = data.frame(row.names = c(1:nc),
								comparison = cmpid, rel.effect = pd, lower = bfn.lower, upper=bfn.upper,
								t.value = t.bf, p.adj = p.bfn, p.raw = p.n), Mult.Distribution = data.frame(Quantile = z.bfn,
								p.Value.global = min(p.bfn)), Correlation = Correlation)
				Asymptotic.Method <- "Multivariate Normal Distribution"
			}, mult.t = {
				x.werte = cbind(lower.bft, pd, upper.bft)
				result <- list(weight.matrix = weight, Data.Info = data.info,
						Analysis.of.relative.effects = data.frame(row.names = c(1:nc),
								comparison = cmpid, rel.effect = pd, lower = bft.lower, upper=bft.upper,
								t.value = t.bf, p.adj = p.bft,p.raw = p.t), Mult.Distribution = data.frame(Quantile = z.bft,
								p.Value.global = min(p.bft), d.f. = df.sw), Correlation = Correlation)
				Asymptotic.Method <- paste("Multi t - Distribution with d.f.= ",
						round(df.sw, 4))
			})
	if (plot.simci == TRUE) {
		test <- matrix(c(1:nc), ncol = nc, nrow = nc)
		angaben <- c(cmpid)
		angaben <- matrix(c(angaben), ncol = nc, nrow = nc)
		k <- c(1:nc)
		plot(x.werte[, 2], k, xlim = c(0, 1), axes = FALSE, type = "p",
				pch = 15, xlab = "", ylab = "")
		abline(v = 0.5, col = "red", lty = 1, lwd = 2)
		axis(1, at = seq(0, 1, 0.1))
		axis(2, at = test, labels = angaben)
		axis(4, at = test, labels = test)
		points(x = x.werte[, 3], y = test[, 1], pch = upper)
		points(x = x.werte[, 1], y = test[, 1], pch = lower)
		for (i in 1:nc) {
			polygon(c(x.werte[i, 1], x.werte[i, 3]), c(i, i))
		}
		box()
		title(main = c(text.output.KI, paste("Type of Contrast:",
								"", type.of.contrast, sep = ""), paste("Method:",
								"", Asymptotic.Method, sep = "")), ylab = "Comparison",
				xlab = paste("lower", lower, "-----", "p", "------",
						upper, "upper"))
	}
	if (info == TRUE) {
		cat("\n", "", "Nonparametric Multiple Comparison Procedure based on relative contrast effects",
				",", "Type of Contrast", ":", type.of.contrast, "\n",
				"NOTE:", "\n", "*-------------------Weight Matrix------------------*",
				"\n", "-", "Weight matrix for choosen contrast based on all-pairs comparisons",
				"\n", "\n", "*-----------Analysis of relative effects-----------*",
				"\n", "-", "Simultaneous Confidence Intervals for relative effects p(i,j)\n      with confidence level",
				conflevel, "\n", "-", "Method", "=", Asymptotic.Method,
				"\n", "-", "p-Values for ", text.output.p, "\n",
				"\n", "*----------------Interpretation--------------------*",
				"\n", "p(a,b)", ">", "1/2", ":", "b tends to be larger than a",
				"\n", "*--------------Mult.Distribution-------------------*",
				"\n", "-", "Equicoordinate Quantile", "\n", "-", "Global p-Value",
				"\n", "*--------------------------------------------------*",
				"\n")
		
	}
	
	return(result)
	
}

weightMatrix <- function(n,type = c("UserDefined","Tukey","AVE","Dunnett", "Sequen",
				"Changepoint", "Marcus",
				"McDermott", "Williams", "UmbrellaWilliams"), base = 1, contrast.matrix=NULL) {
	a  <- length (n)
	
	n.col <- a*(a-1)/2
	tmp <- expand.grid(1:a, 1:a)
	ind <- tmp[[1]] > tmp[[2]]
	vi <- tmp[[2]][ind]
	vj <- tmp[[1]][ind]
	
	
	type <- match.arg(type)
	
	switch (type,
			UserDefined = {
				if (is.null(contrast.matrix)){stop("Choose a contrast or give a contrast matrix by using 'contrast.matrix = matrix'")}
				#if(any(abs(contrast.matrix)>1)){stop("The contrast weights must be between 0 and 1!")}
				if (ncol(contrast.matrix)!=a){stop("The contrast matrix has more or less columns than samples!")}
				nc<-nrow(contrast.matrix)
				
				for (i in 1:nc){
					places_pos<-(contrast.matrix[i,]>0)
					places_neg<-(contrast.matrix[i,]<0)
					sum_pos<-sum(contrast.matrix[i,][places_pos])
					sum_neg<-sum(contrast.matrix[i,][places_neg])
					if (abs(sum_neg)!= sum(sum_pos)) {
						stop(" Wrong contrast matrix!The sum of negative and positive weights must be equal") }
					
					for ( j in 1:a){
						
						if (contrast.matrix[i,j]<0){ contrast.matrix[i,j]<-contrast.matrix[i,j]/abs(sum_neg)}
						if (contrast.matrix[i,j]>0){ contrast.matrix[i,j]<-contrast.matrix[i,j]/sum_pos}
						
					}
				}
				n.row <- nrow(contrast.matrix)
				
				
				ch <- contrast.matrix
				w<-matrix(0,nrow = n.row , ncol=n.col )
				weight.help2<-matrix(0,nrow = n.row , ncol=n.col )
				weight.help1<-matrix(0,nrow = n.row , ncol=n.col )
				for (i in 1:n.row){
					for (j in 1:n.col){
						
						help <- c(rep(i,n.col))
						a <- help[j]
						b <- vi[j]
						d <- vj[j]
						
						w[i,j] <- ch[a,b] *ch[a,d]
						
						if (w[i,j]>0){w[i,j]<-0}
						if (ch[a,b] > 0 && ch[a,d] < 0) {weight.help2[i,j]<-1}
						if(ch[a,b] < 0 && ch[a,d] > 0){weight.help1[i,j] <-1}
					}
				}
				w1<-w*weight.help1
				w2<-w*weight.help2
				w<--1*cbind(w1,w2)
				w
				
			},
			
			AVE = {
				
				
				n.row <- a
				
				ch <- contrMat(n, type="AVE")
				w<-matrix(0,nrow = n.row , ncol=n.col )
				weight.help2<-matrix(0,nrow = n.row , ncol=n.col )
				weight.help1<-matrix(0,nrow = n.row , ncol=n.col )
				for (i in 1:n.row){
					for (j in 1:n.col){
						
						help <- c(rep(i,n.col))
						a <- help[j]
						b <- vi[j]
						d <- vj[j]
						
						w[i,j] <- ch[a,b] *ch[a,d]
						
						if (w[i,j]>0){w[i,j]<-0}
						if (ch[a,b] > 0 && ch[a,d] < 0) {weight.help2[i,j]<-1}
						if(ch[a,b] < 0 && ch[a,d] > 0){weight.help1[i,j] <-1}
					}
				}
				w1<-w*weight.help1
				w2<-w*weight.help2
				w<--1*cbind(w1,w2)
				w
			},
			
			Changepoint = {
				
				n.row <- a - 1
				
				ch <- contrMat(n, type = "Changepoint")
				w<-matrix(0,nrow = n.row , ncol=n.col )
				weight.help2<-matrix(0,nrow = n.row , ncol=n.col )
				weight.help1<-matrix(0,nrow = n.row , ncol=n.col )
				for (i in 1:n.row){
					for (j in 1:n.col){
						
						help <- c(rep(i,n.col))
						a <- help[j]
						b <- vi[j]
						d <- vj[j]
						
						w[i,j] <- ch[a,b] *ch[a,d]
						
						if (w[i,j]>0){w[i,j]<-0}
						if (ch[a,b] > 0 && ch[a,d] < 0) {weight.help2[i,j]<-1}
						if(ch[a,b] < 0 && ch[a,d] > 0){weight.help1[i,j] <-1}
					}
				}
				w1<-w*weight.help1
				w2<-w*weight.help2
				w<--1*cbind(w1,w2)
				w
			},
			
			Dunnett = {
				
				n.row <- a - 1
				
				ch <- contrMat(n, type = "Dunnett", base = base)
				w<-matrix(0,nrow = n.row , ncol=n.col )
				weight.help2<-matrix(0,nrow = n.row , ncol=n.col )
				weight.help1<-matrix(0,nrow = n.row , ncol=n.col )
				for (i in 1:n.row){
					for (j in 1:n.col){
						
						help <- c(rep(i,n.col))
						a <- help[j]
						b <- vi[j]
						d <- vj[j]
						
						w[i,j] <- ch[a,b] *ch[a,d]
						
						if (w[i,j]>0){w[i,j]<-0}
						if (ch[a,b] > 0 && ch[a,d] < 0) {weight.help2[i,j]<-1}
						if(ch[a,b] < 0 && ch[a,d] > 0){weight.help1[i,j] <-1}
					}
				}
				w1<-w*weight.help1
				w2<-w*weight.help2
				w<--1*cbind(w1,w2)
				w
			},
			
			Sequen =  {
				
				n.row <- a - 1
				
				ch <- contrMat(n, type = "Sequen", base)
				w<-matrix(0,nrow = n.row , ncol=n.col )
				weight.help2<-matrix(0,nrow = n.row , ncol=n.col )
				weight.help1<-matrix(0,nrow = n.row , ncol=n.col )
				for (i in 1:n.row){
					for (j in 1:n.col){
						
						help <- c(rep(i,n.col))
						a <- help[j]
						b <- vi[j]
						d <- vj[j]
						
						w[i,j] <- ch[a,b] *ch[a,d]
						
						if (w[i,j]>0){w[i,j]<-0}
						if (ch[a,b] > 0 && ch[a,d] < 0) {weight.help2[i,j]<-1}
						if(ch[a,b] < 0 && ch[a,d] > 0){weight.help1[i,j] <-1}
					}
				}
				w1<-w*weight.help1
				w2<-w*weight.help2
				w<--1*cbind(w1,w2)
				w
			},
			
			Marcus = {
				
				n.row <- a*(a-1)/2
				
				ch <- contrMat(n, type = "Marcus", base)
				w<-matrix(0,nrow = n.row , ncol=n.col )
				weight.help2<-matrix(0,nrow = n.row , ncol=n.col )
				weight.help1<-matrix(0,nrow = n.row , ncol=n.col )
				for (i in 1:n.row){
					for (j in 1:n.col){
						
						help <- c(rep(i,n.col))
						a <- help[j]
						b <- vi[j]
						d <- vj[j]
						
						w[i,j] <- ch[a,b] *ch[a,d]
						
						if (w[i,j]>0){w[i,j]<-0}
						if (ch[a,b] > 0 && ch[a,d] < 0) {weight.help2[i,j]<-1}
						if(ch[a,b] < 0 && ch[a,d] > 0){weight.help1[i,j] <-1}
					}
				}
				w1<-w*weight.help1
				w2<-w*weight.help2
				w<--1*cbind(w1,w2)
				w
			},
			
			McDermott = {
				
				n.row <- a -1
				
				ch <- contrMat(n, type = "McDermott", base)
				w<-matrix(0,nrow = n.row , ncol=n.col )
				weight.help2<-matrix(0,nrow = n.row , ncol=n.col )
				weight.help1<-matrix(0,nrow = n.row , ncol=n.col )
				for (i in 1:n.row){
					for (j in 1:n.col){
						
						help <- c(rep(i,n.col))
						a <- help[j]
						b <- vi[j]
						d <- vj[j]
						
						w[i,j] <- ch[a,b] *ch[a,d]
						
						if (w[i,j]>0){w[i,j]<-0}
						if (ch[a,b] > 0 && ch[a,d] < 0) {weight.help2[i,j]<-1}
						if(ch[a,b] < 0 && ch[a,d] > 0){weight.help1[i,j] <-1}
					}
				}
				w1<-w*weight.help1
				w2<-w*weight.help2
				w<--1*cbind(w1,w2)
				w
			},
			
			Williams = {
				
				n.row <- a - 1
				
				ch <- contrMat(n, type = "Williams")
				w<-matrix(0,nrow = n.row , ncol=n.col )
				weight.help2<-matrix(0,nrow = n.row , ncol=n.col )
				weight.help1<-matrix(0,nrow = n.row , ncol=n.col )
				for (i in 1:n.row){
					for (j in 1:n.col){
						
						help <- c(rep(i,n.col))
						a <- help[j]
						b <- vi[j]
						d <- vj[j]
						
						w[i,j] <- ch[a,b] *ch[a,d]
						
						if (w[i,j]>0){w[i,j]<-0}
						if (ch[a,b] > 0 && ch[a,d] < 0) {weight.help2[i,j]<-1}
						if(ch[a,b] < 0 && ch[a,d] > 0){weight.help1[i,j] <-1}
					}
				}
				w1<-w*weight.help1
				w2<-w*weight.help2
				w<--1*cbind(w1,w2)
				w
			} ,
			
			Tukey = {
				
				
				n.row <- a*(a-1)/2
				
				ch <- contrMat(n, type = "Tukey", base)
				w<-matrix(0,nrow = n.row , ncol=n.col )
				weight.help2<-matrix(0,nrow = n.row , ncol=n.col )
				weight.help1<-matrix(0,nrow = n.row , ncol=n.col )
				for (i in 1:n.row){
					for (j in 1:n.col){
						
						help <- c(rep(i,n.col))
						a <- help[j]
						b <- vi[j]
						d <- vj[j]
						
						w[i,j] <- ch[a,b] *ch[a,d]
						
						if (w[i,j]>0){w[i,j]<-0}
						if (ch[a,b] > 0 && ch[a,d] < 0) {weight.help2[i,j]<-1}
						if(ch[a,b] < 0 && ch[a,d] > 0){weight.help1[i,j] <-1}
					}
				}
				w1<-w*weight.help1
				w2<-w*weight.help2
				w<--1*cbind(w1,w2)
				w
			},
			UmbrellaWilliams = {
				
				n.row <- a*(a-1)/2
				
				ch <- contrMat(n, type = "UmbrellaWilliams", base)
				w<-matrix(0,nrow = n.row , ncol=n.col )
				weight.help2<-matrix(0,nrow = n.row , ncol=n.col )
				weight.help1<-matrix(0,nrow = n.row , ncol=n.col )
				for (i in 1:n.row){
					for (j in 1:n.col){
						
						help <- c(rep(i,n.col))
						a <- help[j]
						b <- vi[j]
						d <- vj[j]
						
						w[i,j] <- ch[a,b] *ch[a,d]
						
						if (w[i,j]>0){w[i,j]<-0}
						if (ch[a,b] > 0 && ch[a,d] < 0) {weight.help2[i,j]<-1}
						if(ch[a,b] < 0 && ch[a,d] > 0){weight.help1[i,j] <-1}
					}
				}
				w1<-w*weight.help1
				w2<-w*weight.help2
				w<--1*cbind(w1,w2)
				w
			}
	)
	ch.out <- matrix(c(ch),nrow=n.row)
	
	result <- list(weight.matrix = w, weight.help1 = weight.help1,weight.help2 = weight.help2,contrast.matrix = ch.out)
	result
}

nparcomp.wrapper <- function(model, data, hypotheses, alpha, alternative, asy.method) {
	
	control <- NULL
	type <- NULL
	contrastMatrix <- NULL	
	
	if (hypotheses %in% c("Tukey", "Dunnett", "Sequen",
					"Williams", "Changepoint", "AVE", "McDermott", "Marcus", "UmbrellaWilliams")) {
		type <- hypotheses
	} else {
		type <- "UserDefined" 
		contrastMatrix <- hypotheses
	}	
	
	result <- nparcomp(formula=formula(model), 
			data, 
			type = type,
			control = control, 
			conflevel = 1-alpha,
			alternative,
			rounds = 3, 
			correlation = TRUE,
			asy.method,
			plot.simci = FALSE,
			info = TRUE,
			contrastMatrix = contrastMatrix)
	pvalues <- result$Analysis.of.relative.effects$p.adj
	rejected1 <- (pvalues < alpha)
	confi <- cbind(result$Analysis.of.relative.effects$rel.effect, result$Analysis.of.relative.effects$lower, result$Analysis.of.relative.effects$upper)
	rownames(confi)<-result$Analysis.of.relative.effects$comparison
	temp <- cbind(confi,pvalues)
	colnames(temp) <- c("Estimate", "Lower","Upper","pValue")
	print(temp)
	return(list(adjPValues=pvalues,
					rejected=rejected1, confIntervals= confi,
					errorControl = new(Class='ErrorControl',type="FWER",alpha=alpha)))
}

mutoss.nparcomp <- function() { return(new(Class="MutossMethod",
					label="Nonparametric relative contrast effects",
					errorControl="FWER",
					callFunction="nparcomp.wrapper",
					output=c("adjPValues", "rejected","confIntervals", "errorControl"),
					info="<h2>Nonparametric relative contrast effects</h2>
							<p>With this function, it is possible to compute nonparametric simultaneous confidence\
							intervals for relative contrast effects in the unbalanced one way layout. Moreover, it computes\
							adjusted p-values. The simultaneous confidence intervals can be computed using\
							multivariate normal distribution, multivariate t-distribution with a Satterthwaite Approximation\
							of the degree of freedom or using multivariate range preserving transformations with Logit or\
							Probit as transformation function. There is no assumption on the underlying distribution function, only\
							that the data have to be at least ordinal numbers.</p> 
							<p></p>
						  <h3>Reference:</h3>
							<ul>
							  <li>Konietschke, F. \"<i>Simultane Konfidenzintervalle fuer nichtparametrische relative Kontrasteffekte.</i>\" Dissertation, University of Goettingen, 2009.</li>
							  <li>Konietschke, F., Brunner, E., Hothorn, L.A.  \"<i>Simultaneous confidence intervals for nonparametric relative contrast effects.</i>\" Research report at the University of Hannover, 2009.</li>
							</ul>",
					parameters=list(
							data=list(type="data.frame"),
							model=list(type="ANY"),
							hypotheses=list(type="ANY"),
							alpha=list(type="numeric"),
							alternative=list(type="character", label="Alternative", choices=c("two.sided", "less", "greater")),
							asy.method=list(type="character", label="Asymptotic approx. method", choices=c("logit", "probit", "normal", "mult.t"))
					)
			)) }






