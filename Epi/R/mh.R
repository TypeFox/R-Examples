# Mantel-Haenszel estimate and test

mh <- 
function(cases, denom, compare = 1, levels = c(1, 2), by = NULL, 
	cohort = !is.integer(denom), confidence = 0.9)
{
	ndim <- length(dim(cases))
	edgin <- names(dimnames(cases))
	edgen <- paste("Dimension", 1:ndim)
	if (is.null(edgin))
		edges <- edgen
	else
		edges <- ifelse(edgin == "", edgen, edgin)
	if (is.null(edges)) edges <- rep("", ndim)
	if(length(dim(denom)) != ndim) {
		stop("Cases and Pyrs arrays of unequal dimension")
	}
	if(is.numeric(compare)) {
		comp <- as.integer(compare)
		if(comp < 1 || comp > ndim) {
			stop("Illegal argument: compare")
		}
	}
	else {
		comp <- (1:ndim)[edges == compare]
		if(length(comp) != 1) {
			stop("Illegal argument: compare")
		}
	}
	if(!is.null(by)) {
		if (!is.numeric(by)) {
			mtch <- match(by, edges)
			if (any(is.na(mtch))) {
				stop("Illegal argument: by")
			}
			by <- (1:ndim)[mtch]
		}
		if (any(by < 1 | by > ndim | by == comp)) {
			stop("Illegal argument: by")
		}
	}
	gtxt <- vector("character", 3)
	gtxt[1] <- edges[comp]
	gtxt[2] <- dimnames(cases)[[comp]][levels[1]]
	gtxt[3] <- dimnames(cases)[[comp]][levels[2]]
	ctxt <- edges[-c(comp, by)]
	if (length(ctxt) == 0) ctxt <- as.null()
	others <- (1:ndim)[ - comp]
	select <- function(a, el)
	{
		b <- a[el]
		ifelse(is.na(b), 0, b)
	}
	d1 <- apply(cases, others, select, el = levels[1])
	d2 <- apply(cases, others, select, el = levels[2])
	if(length(d1) == 0 || length(d2) == 0) {
		stop("Illegal argument: levels")
	}
	y1 <- apply(denom, others, select, el = levels[1])
	y2 <- apply(denom, others, select, el = levels[2])
	d <- d1 + d2
	y <- y1 + y2
	if (cohort) {
		qt <- ifelse(y>0, (d1 * y2)/y, 0)
		rt <- ifelse(y>0, (d2 * y1)/y, 0)
		ut <- ifelse(y>0, d1 - ((d * y1)/y), 0)
		vt <- ifelse(y>0, (d * y1 * y2)/(y^2), 0)
	}
	else {
		s1 <- d1 + y1
		s2 <- d2 + y2
		t <- s1 + s2
		qt <- ifelse(t>1, (d1 * y2)/t, 0)
		rt <- ifelse(t>1, (d2 * y1)/t, 0)
		ut <- ifelse(t>1, d1 - ((d * s1)/t), 0)
		vt <- ifelse(t>1, (d * y * s1 * s2)/((t - 1) * (t^2)), 0)
	}
	if(!is.null(by)) {
		if(length(by) < ndim - 1) {
			nby <- match(by, others)
			q <- apply(qt, nby, sum)
			r <- apply(rt, nby, sum)
			u <- apply(ut, nby, sum)
			v <- apply(vt, nby, sum)
		}
		else {
			q <- qt
			r <- rt
			u <- ut
			v <- vt
		}
	}
	else {
		q <- sum(qt)
		r <- sum(rt)
		u <- sum(ut)
		v <- sum(vt)
	}
	rr <- q/r
	se <- sqrt(v/(q * r))
	ch <- (u^2)/v
	ef <- exp( - qnorm((1 - confidence)/2) * se)
	if (cohort) 
		ty <- "Rate ratio"
	else
		ty <- "Odds ratio"
	res <- list(groups = gtxt, control = ctxt, type=ty, 
		q=q, r=r, u=u, v=v,
		ratio = rr, se.log.ratio = se, cl.lower = rr/ef, 
		cl.upper = rr * ef, chisq = ch, p.value = 1 - pchisq(
		ch, 1))
	class(res) <- "mh"
	res
}

print.mh <- 
function(m) {
	cat("\n")
	if (!is.null(m$control)) 
		cat("\nMantel-Haenszel comparison for: ")
	else 
		cat("Comparison for: ")
	cat(m$groups[1], " (", m$groups[2], "versus", m$groups[3], ")\n")
	if (!is.null(m$control))
		cat("controlled for:", m$control, "\n")
	cols <- c(m$type, "CL (lower)", "CL (upper)", 
		"Chisq (1 df)", "p-value")
	nr <- length(m$ratio)
	if (is.array(m$ratio)) {
		dnt <- dimnames(m$ratio)
		size <- dim(m$ratio)
		nw <- length(dnt)
	}
	else {
		rn <- names(m$ratio)
		if (length(rn) > 1) 
			dnt <- list(names(m$ratio))
		else
			dnt <- list("")
		size <- nr
		nw <- 1
	}
	dno <- vector("list", nw+1)
	so <- vector("numeric", nw+1)
	dno[[1]] <- dnt[[1]]
	dno[[2]] <- cols
	so[1] <- size[1]
	so[2] <- 5
	if (nw > 1) for (i in 2:nw) {
		dno[[i+1]] <- dnt[[i]]
		so[i+1] <- size[i]
	} 
	s1 <- size[1]
	tab <- cbind(m$ratio, m$cl.lower, m$cl.upper, m$chisq, m$p.value)
#		as.matrix(m$ratio, nrow=s1), 
#		as.matrix(m$cl.lower, nrow=s1), 
#		as.matrix(m$cl.upper, nrow=s1), 
#		as.matrix(m$chisq, nrow=s1), 
#		as.matrix(m$p.value, nrow=s1) )
                     
	print(array(tab, dim=so, dimnames=dno))
	if (nr > 1) {
		Q <- sum(m$q)
		R <- sum(m$r)
		cat("\nOverall Mantel-Haenszel estimate of", m$type, ":", 
			format(Q/R)) 
		h <- sum(((m$q*R-m$r*Q)^2)/m$v)/(Q*R)
		df <- sum(m$v>0)-1
		cat("\nChi-squared test of heterogeneity:", format(h),
			"(",df," df), p =", format(1-pchisq(h, df)), "\n")
	}
	cat("\n")
}

# Power calculations

mh.power <- function(mh, ratio, alpha=0.05) {
	n.se <- log(ratio)/mh$se.log.ratio
	pnorm(n.se - qnorm(1-alpha/2))
}
	
