#' Propensity Score Analysis using Stratification
#' 
#' 
#' @param Y response variable.
#' @param Tr treatment variable.
#' @param strata strata identifier.
#' @param trim allows for a trimmed mean as outcome measure, where trim is from
#'        0 to .5 (.5 implying median).
#' @param minStrata minimum number of treatment or control unitis within a strata 
#'        to include that strata.
#' @export
psa.strata <- function (Y, Tr, strata, trim = 0, minStrata=5) {
	sizes <- reshape2::melt(table(strata, Tr))
	smallStrata <- unique(sizes[sizes$value < minStrata,]$strata)
	if(length(smallStrata) == length(unique(strata))) {
		stop('Not enough strata to continue.')
	}
	if(length(smallStrata) > 0) {
		# TODO: Probably should print a warning
		strata <- as.character(strata)
		rows <- !strata %in% smallStrata
		Tr <- Tr[rows]
		Y <- Y[rows]
		strata <- strata[rows]
	}
	
	statistic = "mean" #TODO: put back as a parameter?
	tr.mean <- function(x) {
		mean(x, trim = trim)
	}
	if(!trim == 0) {
		statistic <- tr.mean
	}
	n <- length(Y)
	nstrat <- dim(table(strata))
	ct.means <- tapply(Y, list(strata, Tr), statistic)
	ncontrol <- as.data.frame(table(strata, Tr))[1:nstrat, 3]
	ntreat <- as.data.frame(table(strata, Tr))[(nstrat + 1):(2 * nstrat), 3]
	summary.strata <- cbind(ncontrol, ntreat, ct.means)
	wts <- rowSums(summary.strata[, 1:2])
	x <- ct.means[, 1]
	y <- ct.means[, 2]
	d <- y - x
	xr <- range(x)
	yr <- range(y)
	min.xy <- min(xr[1], yr[1])
	max.xy <- max(xr[2], yr[2])
	diff.wtd <- sum(d * wts)/n
	o <- order(Tr)
	ord.strata <- strata[o]
	nc <- table(Tr)[1]
	nt <- table(Tr)[2]
	ord.Y <- Y[o]
	var.0 <- tapply(ord.Y[1:nc], ord.strata[1:nc], var)
	ni.0 <- table(ord.strata[1:nc])
	frac.0 <- var.0/ncontrol
	ncp1 <- nc + 1
	ncpnt <- nc + nt
	var.1 <- tapply(ord.Y[ncp1:ncpnt], ord.strata[ncp1:ncpnt], var)
	ni.1 <- table(ord.strata[ncp1:ncpnt])
	frac.1 <- var.1/ntreat
	se.wtd <- ((sum(frac.0) + sum(frac.1))^0.5)/nstrat
	strat.labels <- sort(unique(strata))
	wtt <- wts/max(wts)
	wtss <- wts/n
	mnx <- sum(x * wtss)
	mny <- sum(y * wtss)
	mnxy <- (mnx + mny)/2
	dimnamez <- NULL
	dimnamezz <- unlist(dimnames(ct.means)[2])
	dimnamez <- unlist(dimnames(ct.means)[2])
	Cname <- unlist(dimnamez)[1]
	Tname <- unlist(dimnamez)[2]
	C.wtd <- mnx
	T.wtd <- mny
	approx.t <- diff.wtd/se.wtd
	df <- n - 2 * nstrat
	colnames(summary.strata) <- c(paste("n.", dimnamezz[1], 
			sep = ""), paste("n.", dimnamezz[2], sep = ""), paste("means.", 
			dimnamezz[1], sep = ""), paste("means.", dimnamezz[2], 
			sep = ""))
	out <- list(summary.strata, C.wtd, T.wtd, diff.wtd, se.wtd, 
				approx.t = approx.t, df = df)
	names(out) <- c("summary.strata", paste("wtd.Mn.", Cname, 
					sep = ""), paste("wtd.Mn.", Tname, sep = ""), "ATE", 
					"se.wtd", "approx.t", "df")
	ci.diff <- -diff.wtd
	ci <- c(ci.diff - qt(0.975, df) * se.wtd, ci.diff + qt(0.975, df) * se.wtd)
	ci.out <- -ci[c(2, 1)]
	out <- list(summary.strata, C.wtd, T.wtd, diff.wtd, se.wtd, approx.t, df, ci.out)
	names(out) <- c("summary.strata", paste("wtd.Mn.", Cname, 
					sep = "", collapse = ""), paste("wtd.Mn.", Tname, 
					sep = ""), "ATE", "se.wtd", "approx.t", "df", "CI.95")
	return(out)
}
