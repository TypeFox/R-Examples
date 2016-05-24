.onAttach <- function(lib, pkg) {
	packageStartupMessage(sprintf("Package %s (%s) loaded.
To cite, see citation(\"%s\")\n", pkg, packageDescription(pkg)$Version, pkg))
}

setClass(Class = "mederrData", representation(data = "data.frame", size = "numeric", numi = "numeric", numj = "numeric"))

setMethod("initialize", "mederrData",
		function(.Object, data = data.frame(NA)) {
			tmp <- as.data.frame(data)
			.Object@data <- tmp
			names(.Object@data) <- c("y", "N", "groupi", "groupj")
			.Object@size <- dim(tmp)[1]
			.Object@numi <- nlevels(as.factor(tmp[, 3]))
			.Object@numj <- nlevels(as.factor(tmp[, 4]))
			.Object
		}
)

setClass(Class = "mederrFit", representation(thetai = "matrix", deltaj = "matrix", gamma = "numeric", sigma2 = "numeric", tau2 = "numeric", 
				p.acc.i = "numeric", p.acc.j = "numeric", tune.theta = "numeric", tune.delta = "numeric", k = "numeric", eta = "numeric"))

setMethod("initialize", "mederrFit",
		function(.Object, thetai = matrix(NA), deltaj = matrix(NA), gamma = numeric(0), sigma2 = numeric(0), tau2 = numeric(0), p.acc.i = numeric(0), 
				p.acc.j = numeric(0), tune.theta = numeric(0), tune.delta = numeric(0), k = numeric(0), eta = numeric(0)) {
			.Object@thetai <- thetai
			.Object@deltaj <- deltaj
			.Object@gamma <- gamma
			.Object@sigma2 <- sigma2
			.Object@tau2 <- tau2
			.Object@p.acc.i <- p.acc.i 
			.Object@p.acc.j <- p.acc.j
			.Object@tune.theta <- tune.theta
			.Object@tune.delta <- tune.delta
			.Object@k <- k
			.Object@eta <- eta
			.Object
		}
)

setClass(Class = "mederrResample", representation(log.ir = "array", samp = "array", A = "array", t.new = "array", t.old = "numeric", grd = "list"))

setMethod("initialize", "mederrResample",
		function(.Object, log.ir = array(NA), samp = array(NA), A = array(NA), t.new = array(NA), t.old = numeric(0), grd = list()) {
			.Object@log.ir <- log.ir
			.Object@samp <- samp
			.Object@A <- A
			.Object@t.new <- t.new
			.Object@t.old <- t.old
			.Object@grd <- grd
			.Object
		}
)

setMethod("summary",
		"mederrData",
		function(object, plot = FALSE, ...) {
			if (plot) {
				plot(object)
			}
			y <- object@data
			err.dim <- object@numi
			harm.dim <- 2
			hosp.dim <- object@numj
			
			y.sort <- y[sort(y$groupj, index.return = TRUE)$ix, ]
			tbl <- array(NA, c(harm.dim, err.dim, hosp.dim))
			for (j in 1:hosp.dim) tbl[, , j] <- t(y.sort[(1 + err.dim*(j - 1)):(err.dim*j), 1:2])
			
			marg.tmp.N <- tbl[2, , ]
			marg.hosp.N <- colSums(marg.tmp.N)
			marg.err.N <- rowSums(marg.tmp.N)
			marg.tmp.y <- tbl[1, , ]
			marg.hosp.y <- colSums(marg.tmp.y)
			marg.err.y <- rowSums(marg.tmp.y)
			marg.hosp.p <- marg.hosp.y/marg.hosp.N
			marg.err.p <- marg.err.y/marg.err.N

			err <- list("Error profiles - # of reports" = summary(marg.err.N), "Error profiles - # of reports with harm" = summary(marg.err.y),
					"Error profiles - raw rates of harm" = summary(marg.err.p))
			hosp <- list("Hospitals - # of reports" = summary(marg.hosp.N), "Hospitals - # of reports with harm" = summary(marg.hosp.y),
					"Hospitals - raw rates of harm" = summary(marg.hosp.p))
			print(paste("Number of error profiles = ", err.dim, sep = ""), quote = FALSE)
			print(err, quote = FALSE)
			print(paste("Number of hospitals = ", hosp.dim, sep = ""), quote = FALSE)
			print(hosp, quote = FALSE)
			print(paste("Total number of reports = ", sum(marg.hosp.N), sep = ""), quote = FALSE)
		}
)

setMethod("plot",
		signature(x = "mederrData", y = "missing"),
		function(x, nbins.err = 5, nbins.hosp = 5, jittering = .05, ...) {
			y <- x@data
			err.dim <- x@numi
			harm.dim <- 2
			hosp.dim <- x@numj
			
			y.sort <- y[sort(y$groupj, index.return = TRUE)$ix, ]
			tbl <- array(NA, c(harm.dim, err.dim, hosp.dim))
			for (j in 1:hosp.dim) tbl[, , j] <- t(y.sort[(1 + err.dim*(j - 1)):(err.dim*j), 1:2])
			
			marg.tmp.N <- tbl[2, , ]
			marg.hosp.N <- colSums(marg.tmp.N)
			marg.err.N <- rowSums(marg.tmp.N)
			marg.tmp.y <- tbl[1, , ]
			marg.hosp.y <- colSums(marg.tmp.y)
			marg.err.y <- rowSums(marg.tmp.y)
			marg.hosp.p <- marg.hosp.y/marg.hosp.N
			marg.err.p <- marg.err.y/marg.err.N
			
			if (length(nbins.err) == 1) {
				nbins <- nbins.err + 1
				nbins.e <- nbins
				if (nbins >= max(marg.err.N)/2) {
					stop("Too many bins for the available error profile data. Try with a smaller 'nbins.err' argument value.")
				}
				bins.err <- round(seq(min(marg.err.N), max(marg.err.N), length.out = nbins))
			}
			else if (length(nbins.err) > 1) {
				nbins <- length(nbins.err)
				nbins.e <- nbins
				bins.err <- nbins.err
				if ((bins.err[1] > min(marg.err.N)) | (bins.err[nbins] < max(marg.err.N))) {
					ANSWER <- readline("Provided 'nbins.err' do not cover the observed number of error profiles. Continue plotting anyway? ")
					if ((substr(ANSWER, 1, 1) == "N") | (substr(ANSWER, 1, 1) == "n")) {
						stop("Execution stopped!")
					}
				}
			}
			bins.err.lbl <- character(length(bins.err) - 1)
			for (i in 2:length(bins.err)) {
				if (i == 2) {
					bins.err.lbl[i - 1] <- paste("[", bins.err[i - 1], ",", bins.err[i], "]", sep = "")
				}
				else {
					bins.err.lbl[i - 1] <- paste("(", bins.err[i - 1], ",", bins.err[i], "]", sep = "")
				}
			}
			bins.errors <- numeric(err.dim)
			tmp.err <- data.frame(N = marg.err.N, idx = 1:length(marg.err.N))
			tmp.err <- tmp.err[sort(tmp.err$N, index.return = TRUE)$ix, ]
			j <- 1
			for (i in 2:length(bins.err)) {
				while ((tmp.err[j, 1] <= bins.err[i]) & (j <= err.dim)) {
					bins.errors[j] <- (i - 1)
					j <- j + 1
				}
			}
			tmp.err <- cbind(tmp.err, nbins = bins.errors)
			tmp.err <- tmp.err[sort(tmp.err$idx, index.return = TRUE)$ix, ]
			err <- data.frame(y = marg.err.y, N = marg.err.N, p = marg.err.p, nbins = tmp.err$nbins)
			
			if (length(nbins.hosp) == 1) {
				nbins <- nbins.hosp + 1
				nbins.h <- nbins
				if (nbins >= max(marg.hosp.N)/2) {
					stop("Too many bins for the available hospital data. Try with a smaller 'nbins.hosp' argument value.")
				}
				bins.hosp <- round(seq(min(marg.hosp.N), max(marg.hosp.N), length.out = nbins))
			}
			else if (length(nbins.hosp) > 1) {
				nbins <- length(nbins.hosp)
				nbins.h <- nbins
				bins.hosp <- nbins.hosp
				if ((bins.hosp[1] > min(marg.hosp.N)) | (bins.hosp[nbins] < max(marg.hosp.N))) {
					ANSWER <- readline("Provided 'nbins.hosp' do not cover the observed number of hospitals reports. Continue plotting anyway? ")
					if ((substr(ANSWER, 1, 1) == "N") | (substr(ANSWER, 1, 1) == "n")) {
						stop("Execution stopped!")
					}
				}
			}
			bins.hosp.lbl <- character(length(bins.hosp) - 1)
			for (i in 2:length(bins.hosp)) {
				if (i == 2) {
					bins.hosp.lbl[i - 1] <- paste("[", bins.hosp[i - 1], ",", bins.hosp[i], "]", sep = "")
				}
				else {
					bins.hosp.lbl[i - 1] <- paste("(", bins.hosp[i - 1], ",", bins.hosp[i], "]", sep = "")
				}
			}
			bins.hospitals <- numeric(hosp.dim)
			tmp.hosp <- data.frame(N = marg.hosp.N, idx = 1:length(marg.hosp.N))
			tmp.hosp <- tmp.hosp[sort(tmp.hosp$N, index.return = TRUE)$ix, ]
			j <- 1
			for (i in 2:length(bins.hosp)) {
				while ((tmp.hosp[j, 1] <= bins.hosp[i]) & (j <= hosp.dim)) {
					bins.hospitals[j] <- (i - 1)
					j <- j + 1
				}
			}
			tmp.hosp <- cbind(tmp.hosp, nbins = bins.hospitals)
			tmp.hosp <- tmp.hosp[sort(tmp.hosp$idx, index.return = TRUE)$ix, ]
			hosp <- data.frame(y = marg.hosp.y, N = marg.hosp.N, p = marg.hosp.p, nbins = tmp.hosp$nbins)
			
			par(mfrow = c(2, 2))
			
			boxplot(jitter(p, amount = jittering) ~ nbins, data = err, subset = (nbins > 0), main = "Error profiles", 
					names = bins.err.lbl[as.numeric(names(table(err$nbins)))], ylab = "p", cex.axis = .7, xaxt = "n")
			axis(1, at = seq(1, nbins.e, by = 2), labels = (bins.err.lbl[as.numeric(names(table(err$nbins)))])[seq(1, nbins.e, by = 2)], cex.axis = .7)
			boxplot(jitter(p, amount = jittering) ~ nbins, data = hosp, subset = (nbins > 0), main = "Hospitals",
					names = bins.hosp.lbl[as.numeric(names(table(hosp$nbins)))], ylab = "p", cex.axis = .7, xaxt = "n")
			axis(1, at = seq(1, nbins.h, by = 2), labels = (bins.hosp.lbl[as.numeric(names(table(hosp$nbins)))])[seq(1, nbins.h, by = 2)], cex.axis = .7)
			
			barplot(height = as.numeric(table(err$nbins[err$nbins > 0])), space = 0, ylab = "Frequency", xlab = expression(N['i+']))
			axis(1, at = seq(.5, nbins.e, by = 2), labels = (bins.err.lbl[as.numeric(names(table(err$nbins)))])[seq(1, nbins.e, by = 2)], cex.axis = .7)
			barplot(height = as.numeric(table(hosp$nbins[hosp$nbins > 0])), space = 0, ylab = "Frequency", xlab = expression(N['+j']))
			axis(1, at = seq(.5, nbins.h, by = 2), labels = (bins.hosp.lbl[as.numeric(names(table(hosp$nbins)))])[seq(1, nbins.h, by = 2)], cex.axis = .7)
		}
)

setMethod("summary",
		"mederrFit",
		function(object, ebdm.fit, dat, K = 10, ...) {
			if (nargs() < 3)  {
				stop("The first three arguments should be an object of class 'mederrFit', an estimated EBDM model and an object of class 'mederrData'. Type '?mederrFit' for further help.")
			}
			ni <- dat@numi

			r.bhm <- bayes.rank(object)
			r.ebdm <- rank(EBGM(ebdm.fit)[1:ni])

			oi <- order(r.bhm, decreasing = TRUE)
			log.odds <- object@gamma + object@thetai
			p <- colMeans(exp(log.odds)/(1 + exp(log.odds)))
			w <- oi[1:K]

			out <- data.frame(round(r.bhm[w], 1), r.ebdm[w], (100*round(p, 4)[w]), dat@data$y[w], dat@data$N[w], dat@data$groupi[w], dat@data$groupj[w])
			names(out) <- c("BHM", "EBDM", "100*p_hat", "y", "N", "Error profile", "Hospital")
			out
		}
)

setMethod("plot",
		signature(x = "mederrFit", y = "mederrFit"),
		function(x, y, dat, ...) {
			fit <- x
			fit2 <- y
			
			lo.old <- colMeans(fit@gamma + fit@thetai)
			lo.new <- colMeans(fit2@gamma + fit2@thetai)
			
			ni <- dat@numi
			
			par(mfrow = c(2, 1))
			
			par(mar = c(4, 4, 2, 1))
			plot(c(lo.new, lo.old), rep(1:2, each = ni), ylim = c(.25, 2.75), xlab = "Log odds of harm", ylab = "", main = "Error profiles", yaxt = "n")
			segments(lo.new, rep(1, ni), lo.old, rep(2, ni))
			axis(2, at = 1:2, labels = c("", ""))
			mtext(bquote(k == infinity), side = 2, line = 1, at = 2.15, las = 1)
			mtext(bquote(eta == 1), side = 2, line = 1.1, at = 1.85, las = 1)
			if (y@k == Inf) {
				mtext(bquote(k == infinity), side = 2, line = 1.3, at = 1.15, las = 1)
			}
			else {
				mtext(bquote(k == .(y@k)), side = 2, line = 1.3, at = 1.15, las = 1)
			}
			mtext(bquote(eta == .(y@eta)), side = 2, line = 1, at = 0.85, las = 1)
			abline(h = 1)
			abline(h = 2)
			
			r.old <- bayes.rank(x)
			r.new <- bayes.rank(y)
			
			par(mar = c(4,4,2,1))
			plot(c(r.new, r.old), rep(1:2, each = ni), ylim = c(.25, 2.75), xlab = "Optimal Bayesian rank", ylab = "", main = "Error profiles", yaxt = "n")
			segments(r.new, rep(1, ni), r.old, rep(2, ni))
			axis(2, at = 1:2, labels = c("",""))
			mtext(bquote(k == infinity), side = 2, line = 1, at = 2.15, las = 1)
			mtext(bquote(eta == 1), side = 2, line = 1.1, at = 1.85, las = 1)
			if (y@k == Inf) {
				mtext(bquote(k == infinity), side = 2, line = 1.3, at = 1.15, las = 1)
			}
			else {
				mtext(bquote(k == .(y@k)), side = 2, line = 1.3, at = 1.15, las = 1)
			}
			mtext(bquote(eta == .(y@eta)), side = 2, line = 1, at = 0.85, las = 1)
			abline(h = 1)
			abline(h = 2)
		}
)

setMethod("plot",
		signature(x = "mederrResample", y = "missing"),
		function(x, ...) {
			lp <- logunpost(x)
			interaction.plot(col(lp), row(lp), lp, xaxt = "n", xlab = bquote(eta), lty = 1:length(x@grd$k), legend = FALSE,
					ylab = "Log unnormalized marginal density")
			axis(1, at = 1:length(x@grd$eta), labels = x@grd$eta)
			legend((length(x@grd$eta) - 1), max(lp), legend = tapply(x@grd$k, 1:length(x@grd$k), function(z) paste("k = ", z, sep = "")),
					lty = (1:length(x@grd$k)))
		}
)
