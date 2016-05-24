"epi.prev" <- function(pos, tested, se, sp, method = "wilson", conf.level = 0.95){

   # Apparent prevalence:
   ap.p <- pos / tested
   if(method == "c-p") ap.cl = .bin.conf(pos, tested, method = "e", alpha = 1 - conf.level)[2:3]
   else if (method == "sterne") ap.cl = .sterne.int(pos, tested, alpha = 1 - conf.level)
   else if (method == "blaker") ap.cl = .blakerci(pos, tested, conf.level)
   else if (method == "wilson") ap.cl = .bin.conf(pos, tested, method = "w", alpha = 1 - conf.level)[2:3]
   else stop('Valid methods are "c-p", "sterne", "blaker", or "wilson"')
   
   # True prevalence:
   if(method == "c-p") tp.cl = .bin.conf(pos, tested, method = "e", alpha = 1 - conf.level)[2:3]
   else if (method == "sterne") tp.cl = .sterne.int(pos, tested, alpha = 1 - conf.level)
   else if (method == "blaker") tp.cl = .blakerci(pos, tested, conf.level)
   else if (method == "wilson") tp.cl = .bin.conf(pos, tested, method = "w", alpha = 1 - conf.level)[2:3]
   else stop('Valid methods are "c-p", "sterne", "blaker", or "wilson"')

   tp.p <- (ap.p + sp - 1) / (se + sp - 1)
   tp.p[tp.p < 0] <- 0
   tp.p[tp.p > 1] <- 1
   adj.cl <- (tp.cl + sp - 1) / (se + sp - 1) 
   adj.cl <- pmax(adj.cl, c(0, 0))
   adj.cl <- pmin(adj.cl, c(1, 1))

   result.01 <- data.frame(est = ap.p, lower = ap.cl[1], upper = ap.cl[2])
   result.02 <- data.frame(est = tp.p, lower = adj.cl[1], upper = adj.cl[2])
     
   rval <- list(ap = result.01, tp = result.02)
   return(rval)
}

# library(Hmisc)

# -----------------------------------
# Blaker's interval (by Helge Blaker). Computes the Blaker exact CI (Canadian J. Stat 2000) for a binomial success probability for x successes out of n trials with confidence coefficient = conf.level.
# uses acceptbin function.

.blakerci <- function(x, n, conf.level, tolerance = 1e-04){
   lower = 0
   upper = 1
   if (x != 0){lower = qbeta((1 - conf.level) / 2, x, n - x + 1)
    while (.acceptbin(x, n, lower + tolerance) < (1 - conf.level))
    lower = lower + tolerance
   }
   if (x != n){upper = qbeta(1 - (1 - conf.level) / 2, x + 1, n - x)
    while (.acceptbin(x, n, upper - tolerance) < (1 - conf.level))
    upper = upper - tolerance
   }
   c(lower, upper)
   }


.acceptbin = function(x, n, p){
   # Computes the Blaker acceptability of p when x is observed and X is bin(n, p)
   p1 = 1 - pbinom(x - 1, n, p)
   p2 = pbinom(x, n, p)
   a1 = p1 + pbinom(qbinom(p1, n, p) - 1, n, p)
   a2 = p2 + 1 - pbinom(qbinom(1 - p2, n, p), n, p)
   return(min(a1,a2))
}

# -----------------------------------
# Exact confidence intervals
# -----------------------------------

.sterne.int <- function(x, n, alpha = 0.05, del = 10^-5){
   logit <- function(p){log(p / (1 - p))}
   invlogit <- function(y){exp(y) / (1 + exp(y))}
   theta <- function(k, x, n){(lchoose(n, x) - lchoose(n, k)) / (k - x)}
   Feta <- function(x, eta){pbinom(x, n, invlogit(eta))}

# The function piXeta(x, eta) automatically accounts for the fact that if k_alpha(X) = min(J) then a_alpha^st(X) = a_alpha(X)
.piXeta <- function(x, eta){
	 if (invlogit(eta) >= 1){f <- 0} else {
	 J <- c(0:(x - 1),(x + 1):n)

	 # on (-infinity, theta_0]
	 t1 <- theta(0, x, n)
   if (is.na(t1) != 1 && eta <= t1){f <- 1 - Feta(x - 1, eta)}

	 # on [theta_0,mode]
	 k1 <- J[J < (x - 1)]

	 if (length(k1) > 0){
	 the1 <- theta(k1, x, n)
	 the2 <- theta(k1 + 1, x, n)
	 pos <- (the1 <= eta) * (eta < the2)
	 if (sum(pos) > 0){f <- 1 - Feta(x - 1, eta) + Feta(max(k1 * pos), eta)}
    }

	 # mode
	 the1 <- theta(x - 1, x, n)
	 the2 <- theta(x + 1, x, n)
	 if (eta >= the1 && eta <= the2){f <- 1}
   }

	 # on [mode,theta_n]
	 k2 <- J[J > (x + 1)]
	 if (length(k2) > 0){
	    the1 <- theta(k2 - 1, x, n)
	    the2 <- theta(k2, x, n)
	    kre <- sum(k2 * (the1 < eta) * (eta <= the2))
   if (kre > 0){
      f <- 1 - Feta(kre - 1, eta) + Feta(x, eta)}
   }

	 # on [theta_n,infty)
	 t2 <- theta(n, x, n)
	 if (is.na(t2) != 1 && eta >= t2){f <- Feta(x, eta)}
	 f}

   # Lower bound a_alpha^st(X)
   if (x ==0 ){pu <- 0} else {
   J <- c(0:(x - 1), (x + 1):n)
   k1 <- min(J)
   pi1 <- .piXeta(x, theta(k1, x, n))

   # Calculation of k_alpha(X)
   if (pi1 >= alpha){kal <- k1} else {
      k <- x-1
	 while (k1 < k - 1){
		  k2 <- floor((k + k1) / 2)
		  pi2 <- .piXeta(x, theta(k2, x, n))
   if (pi2 >= alpha){k <- k2} 
   else {k1 <- k2}
   }
   kal <- k
   }

   # Calculation of a_alpha^st(X):
   b1 <- theta(kal, x, n)
   pi1 <- 1 - Feta(x - 1, b1) + Feta(kal - 1, b1)
   if (pi1 <= alpha){b <- b1} else {
   b <- max(theta(kal - 1, x, n),logit(del))
   pi <- 1 - Feta(x - 1, b) + Feta(kal - 1, b)
	 while (b1 - b > del || pi1 - pi > del){
   b2 <- (b + b1) / 2
   pi2 <- 1 - Feta(x - 1, b2) + Feta(kal - 1, b2)
   if (pi2 > alpha){
      b1 <- b2
			pi1 <- pi2} else {
			b <- b2
			pi <- pi2}}}
   pu <- invlogit(b)}

   # Upper bound b_alpha^st(X):
   if (x == n){po <- 1} else {
   J <- c(0:(x - 1),(x + 1):n)
   k1 <- max(J)
   pi1 <- .piXeta(x, theta(k1, x, n))

   # Calculation of k_alpha(X):
   if (pi1 >= alpha){kau <- k1} else {
   k <- x + 1
   pi <- 1
   while (k1 > k + 1){
   k2 <- floor((k + k1) / 2)
   pi2 <- .piXeta(x, theta(k2, x, n))
   if (pi2 >= alpha){k <- k2} 
   else {k1 <- k2}
	 }
   kau <- k
   }

   # Calculation of b_alpha^st(X):
   b1 <- theta(kau, x, n)
   pi1 <- 1 - Feta(kau, b1) + Feta(x, b1)

   if (pi1 <= alpha){
   b <- b1
	 po <- pi1} else {
   b <- min(theta(kau + 1, x, n), b1 + n)
   pi <- 1 - Feta(kau, b) + Feta(x, b)
   while (b - b1 > del || pi1 - pi > del){
   b2 <- (b + b1) / 2
   pi2 <- 1 - Feta(kau, b2) + Feta(x, b2)
   if (pi2 > alpha){
   b1 <- b2
	 pi1 <- pi2} else {
   b <- b2
   pi <- pi2}}}
   po <- invlogit(b)}

   # c("a_alpha^St" = pu, "b_alpha^St" = po)
   c(pu, po)
}

.bin.conf <- function (x, n, alpha = 0.05, method = c("wilson", "exact", "asymptotic", "all"), include.x = FALSE, include.n = FALSE, return.df = FALSE){
    method <- match.arg(method)
    bc <- function(x, n, alpha, method) {
        nu1 <- 2 * (n - x + 1)
        nu2 <- 2 * x
        ll <- if (x > 0)
            x/(x + qf(1 - alpha/2, nu1, nu2) * (n - x + 1))
        else 0
        nu1p <- nu2 + 2
        nu2p <- nu1 - 2
        pp <- if (x < n)
            qf(1 - alpha/2, nu1p, nu2p)
        else 1
        ul <- ((x + 1) * pp)/(n - x + (x + 1) * pp)
        zcrit <- -qnorm(alpha/2)
        z2 <- zcrit * zcrit
        p <- x/n
        cl <- (p + z2/2/n + c(-1, 1) * zcrit * sqrt((p * (1 -
            p) + z2/4/n)/n))/(1 + z2/n)
        if (x == 1)
            cl[1] <- -log(1 - alpha)/n
        if (x == (n - 1))
            cl[2] <- 1 + log(1 - alpha)/n
        asymp.lcl <- x/n - qnorm(1 - alpha/2) * sqrt(((x/n) *
            (1 - x/n))/n)
        asymp.ucl <- x/n + qnorm(1 - alpha/2) * sqrt(((x/n) *
            (1 - x/n))/n)
        res <- rbind(c(ll, ul), cl, c(asymp.lcl, asymp.ucl))
        res <- cbind(rep(x/n, 3), res)
        switch(method, wilson = res[2, ], exact = res[1, ], asymptotic = res[3,
            ], all = res, res)
    }
    if ((length(x) != length(n)) & length(x) == 1)
        x <- rep(x, length(n))
    if ((length(x) != length(n)) & length(n) == 1)
        n <- rep(n, length(x))
    if ((length(x) > 1 | length(n) > 1) & method == "all") {
        method <- "wilson"
        warning("method = all will not work with vectors ... setting method to wilson")
    }
    if (method == "all" & length(x) == 1 & length(n) == 1) {
        mat <- bc(x, n, alpha, method)
        dimnames(mat) <- list(c("Exact", "Wilson", "Asymptotic"),
            c("PointEst", "Lower", "Upper"))
        if (include.n)
            mat <- cbind(N = n, mat)
        if (include.x)
            mat <- cbind(X = x, mat)
        if (return.df)
            mat <- as.data.frame(mat)
        return(mat)
    }
    mat <- matrix(ncol = 3, nrow = length(x))
    for (i in 1:length(x)) mat[i, ] <- bc(x[i], n[i], alpha = alpha,
        method = method)
    dimnames(mat) <- list(rep("", dim(mat)[1]), c("PointEst",
        "Lower", "Upper"))
    if (include.n)
        mat <- cbind(N = n, mat)
    if (include.x)
        mat <- cbind(X = x, mat)
    if (return.df)
        mat <- as.data.frame(mat, row.names = NULL)
    mat
}
