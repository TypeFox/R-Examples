##' Calculate cell central moments
##' 
##' Popular central moments include 2 (variance) and 4 (kurtosis).
##' 
##' @param spec list of item models
##' @param params data frame of item parameters, 1 per row
##' @param scores model derived person scores
##' @param m which moment
##' @return moment matrix
##' @docType methods
##' @export
rpf.1dim.moment <- function(spec, params, scores, m) {
  out <- array(dim=c(length(scores), length(spec)))
  for (ix in 1:length(spec)) {
    i <- spec[[ix]]
    prob <- t(rpf.prob(i, params[,ix], scores))  # remove t() TODO
    Escore <- apply(prob, 1, function(r) sum(r * 0:(i@outcomes-1)))
    grid <- t(array(0:(i@outcomes-1), dim=c(i@outcomes, length(scores))))
    out[,ix] <- apply((grid - Escore)^m * prob, 1, sum)
  }
  out
}

##' Calculate residuals
##'
##' @param spec list of item models
##' @param params data frame of item parameters, 1 per row
##' @param responses persons in rows and items in columns
##' @param scores model derived person scores
##' @return residuals
##' @docType methods
##' @export
rpf.1dim.residual <- function(spec, params, responses, scores) {
  Zscore <- array(dim=c(length(scores), length(spec)))
  for (ix in 1:length(spec)) {
    i <- spec[[ix]]
    prob <- t(rpf.prob(i, params[,ix], scores))  # remove t() TODO
    Escore <- apply(prob, 1, function(r) sum(r * 0:(i@outcomes-1)))
    data <- responses[,ix]
    if (is.ordered(data)) {
      data <- unclass(data) - 1
    } else if (is.factor(data)) {
      stop(paste("Column",ix,"is an unordered factor"))
    }
    if (length(data) != length(Escore)) stop("Length mismatch")
    Zscore[,ix] <- data - Escore
  }
  Zscore
}

##' Calculate standardized residuals
##'
##' @param spec list of item models
##' @param params data frame of item parameters, 1 per row
##' @param responses persons in rows and items in columns
##' @param scores model derived person scores
##' @return standardized residuals
##' @docType methods
##' @export
rpf.1dim.stdresidual <- function(spec, params, responses, scores) {
  res <- rpf.1dim.residual(spec, params, responses, scores)
  variance <- rpf.1dim.moment(spec, params, scores, 2)
  res / sqrt(variance)
}

##' Calculate item and person Rasch fit statistics
##'
##' Note: These statistics are only appropriate if all discrimination
##' parameters are fixed equal and items are conditionally independent
##' (see \code{\link{ChenThissen1997}}).  A best effort is made to
##' cope with missing data.
##'
##' Exact distributional properties of these statistics are unknown
##' (Masters & Wright, 1997, p. 112).  For details on the calculation,
##' refer to Wright & Masters (1982, p. 100).
##'
##' The Wilson-Hilferty transformation is biased for less than 25 items.
##' Consider wh.exact=FALSE for less than 25 items.
##'
##' @param spec list of item models
##' @param params matrix of item parameters, 1 per column
##' @param responses persons in rows and items in columns
##' @param scores model derived person scores
##' @param margin for people 1, for items 2
##' @param wh.exact whether to use the exact Wilson-Hilferty transformation
##' @param group spec, params, data, and scores can be provided in a list instead of as arguments
##' @references Masters, G. N. & Wright, B. D. (1997). The Partial
##' Credit Model. In W. van der Linden & R. K. Kambleton (Eds.),
##' \emph{Handbook of modern item response theory}
##' (pp. 101-121). Springer.
##' 
##' Wilson, E. B., & Hilferty, M. M. (1931). The distribution of
##' chi-square. \emph{Proceedings of the National Academy of Sciences of the
##' United States of America,} 17, 684-688.
##' 
##' Wright, B. D. & Masters, G. N. (1982). \emph{Rating Scale
##' Analysis.} Chicago: Mesa Press.
##' @export
##' @examples
##' data(kct)
##' responses <- kct.people[,paste("V",2:19, sep="")]
##' rownames(responses) <- kct.people$NAME
##' colnames(responses) <- kct.items$NAME
##' scores <- kct.people$MEASURE
##' params <- cbind(1, kct.items$MEASURE, logit(0), logit(1))
##' rownames(params) <- kct.items$NAME
##' items<-list()
##' items[1:18] <- rpf.drm()
##' params[,2] <- -params[,2]
##' rpf.1dim.fit(items, t(params), responses, scores, 2, wh.exact=TRUE)
rpf.1dim.fit <- function(spec, params, responses, scores, margin, group=NULL, wh.exact=TRUE) {
    if (!missing(group)) {
        spec <- group$spec
        params <- group$param
        responses <- group$data
        scores <- group$score[,1]  # should not assume first score TODO
    }

  if (any(is.na(responses))) warning("Rasch fit statistics should not be used with missing data")  # true? TODO

  if (dim(params)[2] < 25 && wh.exact) {
    if (missing(wh.exact)) {
      wh.exact <- FALSE
      warning("Consider wh.exact=FALSE for less than 25 items")
    }
  }
  if (dim(params)[2] > 25 && !wh.exact) {
    if (missing(wh.exact)) {
      wh.exact <- TRUE
      warning("Consider wh.exact=TRUE for more than 25 items")
    }
  }

  exclude.col <- c()
  outcomes <- sapply(spec, function(s) s@outcomes)
  for (ix in 1:dim(responses)[2]) {
    kat <- sum(table(responses[,ix]) > 0)
    if (kat != outcomes[ix]) {
      exclude.col <- c(exclude.col, ix)
      warning(paste("Excluding item", colnames(responses)[ix], "because outcomes !=", outcomes[ix]))
    }
  }

  if (length(exclude.col)) {
    responses <- responses[,-exclude.col]
    spec <- spec[-exclude.col]
    params <- params[,-exclude.col]
    outcomes <- outcomes[-exclude.col]
  }

  exclude.row <- c()
  for (ix in 1:dim(responses)[1]) {
    r1 <- sapply(responses[ix,], unclass)
    if (any(is.na(r1))) next
    if (all(r1 == 1) || all(r1 == outcomes)) {
      exclude.row <- c(exclude.row, ix)
      warning(paste("Excluding response", rownames(responses)[ix], "because it is a minimum or maximum"))
    }
  }

  if (length(exclude.row)) {
    responses <- responses[-exclude.row,]
    scores <- scores[-exclude.row]
  }

  na.rm=TRUE
  r.z <- rpf.1dim.stdresidual(spec, params, responses, scores)
  r.var <- rpf.1dim.moment(spec, params, scores,2)
  r.var[is.na(r.z)] <- NA
  r.k <- rpf.1dim.moment(spec, params, scores,4)
  r.k[is.na(r.z)] <- NA

  outfit.var <- r.var
  outfit.var[r.var^2 < 1e-5] <- sqrt(1e-5)
  outfit.n <- apply(r.var, margin, function(l) sum(!is.na(l)))
  outfit.sd <- sqrt(apply(r.k / outfit.var^2, margin, sum, na.rm=na.rm) /
                    outfit.n^2 - 1/outfit.n)
  outfit.sd[outfit.sd > 1.4142] <- 1.4142
  outfit.fudge <- outfit.sd/3

  infit.sd <- sqrt(apply(r.k - r.var^2, margin, sum, na.rm=na.rm)/
    apply(r.var, margin, sum, na.rm=na.rm)^2)
  infit.sd[infit.sd > 1.4142] <- 1.4142
  infit.fudge <- infit.sd/3
  if (!wh.exact) {
    infit.fudge <- 0
    outfit.fudge <- 0
  }

  outfit <- apply(r.z^2, margin, sum, na.rm=na.rm)/
                       apply(r.z, margin, function (l) sum(!is.na(l)))
  outfit.z <- (outfit^(1/3) - 1)*(3/outfit.sd) + outfit.fudge

  infit <- apply(r.z^2 * r.var, margin, sum, na.rm=na.rm)/
                     apply(r.var, margin, sum, na.rm=na.rm)
  infit.z <- (infit^(1/3) - 1)*(3/infit.sd) + infit.fudge

  df <- data.frame(n=outfit.n, infit, infit.z, outfit, outfit.z)
  if (margin == 2) {
    df$name <- colnames(params)
  } else {
    df$name <- rownames(responses)
  }
  df
}

##' Find the point where an item provides mean maximum information
##'
##' WARNING: This function is experimental and may disappear.
##'
##' @param spec an item spec
##' @param iparam an item parameter vector
##' @param grain the step size for numerical integration (optional)
rpf.mean.info1 <- function(spec, iparam, grain=.1) {
  range <- 9
  dim <- spec@factors
  if (dim != 1) stop("Not implemented")
  grid <- seq(-range, range, grain)
  info <- rpf.info(spec, iparam, grid)
  sum(info * grid) / sum(info)
}

##' Find the point where an item provides mean maximum information
##'
##' This is a point estimate of the mean difficulty of items that do
##' not offer easily interpretable parameters such as the Generalized
##' PCM. Since the information curve may not be unimodal, this
##' function integrates across the latent space.
##' 
##' WARNING: This function is experimental and may disappear.
##' 
##' @param spec list of item specs
##' @param param list or matrix of item parameters
##' @param grain the step size for numerical integration (optional)
rpf.mean.info <- function(spec, param, grain=.1) {
  ret <- list()
  for (ix in 1:length(spec)) {
    iparam <- c()
    if (is.list(param)) {
      iparam <- param[[ix]]
    } else {
      iparam <- param[ix,]
    }
    ret[[ix]] <- rpf.mean.info1(spec[[ix]], iparam, grain)
  }
  ret
}

# copied from mirt
collapseCells <- function(On, En, mincell = 1){
		drop <- which(rowSums(is.na(En)) > 0)
		En[is.na(En)] <- 0
					#collapse known upper and lower sparce cells
		if(length(drop) > 0L){
			up <- drop[1L]:drop[length(drop)/2]
			low <- drop[length(drop)/2 + 1L]:drop[length(drop)]
			En[max(up)+1, ] <- colSums(En[c(up, max(up)+1), , drop = FALSE])
			On[max(up)+1, ] <- colSums(On[c(up, max(up)+1), , drop = FALSE])
			En[min(low)-1, ] <- colSums(En[c(low, min(low)-1), , drop = FALSE])
			On[min(low)-1, ] <- colSums(On[c(low, min(low)-1), , drop = FALSE])
			En[c(up, low), ] <- On[c(up, low), ] <- NA
			En <- na.omit(En)
			On <- na.omit(On)
		}
					#collapse accross
		if(ncol(En) > 2L){
			for(j in 1L:(ncol(On)-1L)){
				L <- En < mincell
				sel <- L[,j]
				if(!any(sel)) next
				On[sel, j+1L]  <- On[sel, j] + On[sel, j+1L]
				En[sel, j+1L]  <- En[sel, j] + En[sel, j+1L]
				On[sel, j] <- En[sel, j] <- NA
			}
			sel <- L[,j+1L]
			sel[rowSums(is.na(En[, 1L:j])) == (ncol(En)-1L)] <- FALSE
			put <- apply(En[sel, 1L:j, drop=FALSE], 1, function(x) max(which(!is.na(x))))
			put2 <- which(sel)
			for(k in 1L:length(put)){
				En[put2[k], put[k]] <- En[put2[k], put[k]] + En[put2[k], j+1L]
				En[put2[k], j+1L] <- On[put2[k], j+1L] <- NA
			}
		}
		L <- En < mincell
		L[is.na(L)] <- FALSE
		while(any(L)){
			drop <- c()
			for(j in 1L:(nrow(On)-1L)){
				if(any(L[j,])) {
					On[j+1L, L[j,]] <- On[j+1L, L[j,]] + On[j, L[j,]]
					En[j+1L, L[j,]] <- En[j+1L, L[j,]] + En[j, L[j,]]
					drop <- c(drop, j)
					break
				}
			}
			for(j in nrow(On):2L){
				if(any(L[j,])) {
					On[j-1L, L[j,]] <- On[j-1L, L[j,]] + On[j, L[j,]]
					En[j-1L, L[j,]] <- En[j-1L, L[j,]] + En[j, L[j,]]
					drop <- c(drop, j)
					break
				}
			}
			if(nrow(On) > 4L){
				for(j in 2L:(nrow(On)-1L)){
					if(any(L[j,])){
						On[j+1L, L[j,]] <- On[j+1L, L[j,]] + On[j, L[j,]]
						En[j+1L, L[j,]] <- En[j+1L, L[j,]] + En[j, L[j,]]
						drop <- c(drop, j)
						break
					}
				}
			}
					#drop
			if(!is.null(drop)){
				En <- En[-drop, ]
				On <- On[-drop, ]
			}
			L <- En < mincell
			L[is.na(L)] <- FALSE
		}
        return(list(O=On, E=En))
}

SitemFit1Internal <- function(out) {
	observed <- out$orig.observed
	expected <- out$orig.expected

    mask <- apply(observed, 1, sum) != 0
    observed = observed[mask,,drop=FALSE]
    expected = expected[mask,,drop=FALSE]
	if (!length(observed)) {
		out$statistic <- NA
		out$pval <- NA
		return(out)
	}

	method <- out$method
    if (method == "pearson") {
	if (out$alt) {
		kc <- collapseCells(observed, expected)
	} else {
		kc <- .Call(collapse_wrapper, observed, expected)
	}
        out$observed <- observed <- kc$O
        out$expected <- expected <- kc$E
        mask <- !is.na(expected) & expected!=0
        out$statistic <- sum((observed[mask] - expected[mask])^2 / expected[mask])
	    out$df <- nrow(observed) * (ncol(observed) - 1)
        out$df <- out$df - out$free - sum(is.na(expected));
	out$df[out$df < 1] <- 1
        out$pval <- pchisq(out$statistic, out$df, lower.tail=FALSE, log.p=log)
    } else if (method == "rms") {
      pval <- crosstabTest(observed, expected)
        out$observed <- observed
        out$expected <- expected
        if (log) {
            out$pval <- log(pval)
        } else {
            out$pval <- pval
        }
    } else {
        stop(paste("Method", method, "not recognized"))
    }
    out
}

##' Compute the S fit statistic for 1 item
##'
##' Implements the Kang & Chen (2007) polytomous extension to
##' S statistic of Orlando & Thissen (2000). Rows with
##' missing data are ignored, but see the \code{omit} option.
##'
##' This statistic is good at finding a small number of misfitting
##' items among a large number of well fitting items. However, be
##' aware that misfitting items can cause other items to misfit.
##'
##' Observed tables cannot be computed when data is
##' missing. Therefore, you can optionally omit items with the
##' greatest number of responses missing relative to the item of
##' interest.
##' 
##' Pearson is slightly more powerful than RMS in most cases I
##' examined.
##'
##' Setting \code{alt} to \code{TRUE} causes the tables to match
##' published articles. However, the default setting of \code{FALSE}
##' probably provides slightly more power when there are less than 10
##' items.
##'
##' The name of the test, "S", probably stands for sum-score.
##'
##' @param grp a list with spec, param, mean, cov, and data
##' @param item the item of interest
##' @param free the number of free parameters involved in estimating the item (to adjust the df)
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param method whether to use a pearson or rms test
##' @param log whether to return pvalues in log units
##' @param qwidth the positive width of the quadrature in Z units
##' @param qpoints the number of quadrature points
##' @param alt whether to include the item of interest in the denominator
##' @param omit number of items to omit or a character vector with the names of the items to omit when calculating the observed and expected sum-score tables
##' @param .twotier whether to enable the two-tier optimization
##' @references Kang, T. and Chen, T. T. (2007). An investigation of
##' the performance of the generalized S-Chisq item-fit index for
##' polytomous IRT models. ACT Research Report Series.
##'
##' Orlando, M. and Thissen, D. (2000). Likelihood-Based
##' Item-Fit Indices for Dichotomous Item Response Theory Models.
##' \emph{Applied Psychological Measurement, 24}(1), 50-64.
SitemFit1 <- function(grp, item, free=0, ..., method="pearson", log=TRUE, qwidth=6, qpoints=49L,
		      alt=FALSE, omit=0L, .twotier=TRUE) {
	if (length(list(...)) > 0) {
		stop(paste("Remaining parameters must be passed by name", deparse(list(...))))
	}

	if (missing(qwidth) && !is.null(grp$qwidth)) { qwidth <- grp$qwidth }
	if (missing(qpoints) && !is.null(grp$qpoints)) { qpoints <- grp$qpoints }

    spec <- grp$spec
  c.spec <- lapply(spec, function(m) {
    if (length(m@spec)==0) { stop("Item model",m,"is not implemented") }
    else { m@spec }
  })

    param <- grp$param
    itemIndex <- which(item == colnames(param))

	mask <- rep(TRUE, ncol(param))
	if (!alt) mask[itemIndex] <- FALSE
	omitted <- NULL
	if (is.null(omit)) {
		# OK
	} else if (is.numeric(omit)) {
		omitted <- bestToOmit(grp, omit, item)
	} else if (is.character(omit)) {
		omitted <- omit
	} else {
		stop(paste("Not clear how to interpret omit =", omit))
	}
	mask[match(omitted, colnames(grp$param))] <- FALSE
	iobss <- itemOutcomeBySumScore(grp, mask, itemIndex)
	observed <-iobss$table

  max.param <- max(vapply(spec, rpf.numParam, 0))
  if (nrow(param) < max.param) {
    stop(paste("param matrix must have", max.param ,"rows"))
  }

    Eproportion <- ot2000md(grp, itemIndex, qwidth, qpoints, alt, mask, .twotier)
    if (nrow(Eproportion) != nrow(observed)) {
	    print(Eproportion)
	    stop(paste("Expecting", nrow(observed), "rows in expected matrix"))
    }
    Escale <- matrix(apply(observed, 1, sum), nrow=nrow(Eproportion), ncol=ncol(Eproportion))
    expected <- Eproportion * Escale
	names(dimnames(observed)) <- c("sumScore", "outcome")
	dimnames(expected) <- dimnames(observed)
    out <- list(orig.observed=observed, orig.expected = expected,
		log=log, method=method, n=iobss$n, free=free, alt=alt, omitted=omitted)

	out <- SitemFit1Internal(out)
	out
}

ot2000md <- function(grp, item, width, pts, alt=FALSE, mask, .twotier) {
	if (missing(width)) width <- 6
	if (missing(pts)) pts <- 49L
	.Call(ot2000_wrapper, grp, item, width, pts, alt, mask, .twotier)
}

##' Compute the S fit statistic for a set of items
##'
##' Runs \code{\link{SitemFit1}} for every item and accumulates
##' the results.
##'
##' @param grp a list with spec, param, mean, cov, data, and the free variable pattern
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param method whether to use a pearson or rms test
##' @param log whether to return pvalues in log units
##' @param qwidth the positive width of the quadrature in Z units
##' @param qpoints the number of quadrature points
##' @param alt whether to include the item of interest in the denominator
##' @param omit number of items to omit (a single number) or a list of the length the number of items
##' @param .twotier whether to enable the two-tier optimization
##' @param .parallel whether to take advantage of multiple CPUs (default TRUE)
##' @return
##' a list of output from \code{\link{SitemFit1}}
##' @examples
##' grp <- list(spec=list())
##' grp$spec[1:20] <- rpf.grm()
##' grp$param <- sapply(grp$spec, rpf.rparam)
##' colnames(grp$param) <- paste("i", 1:20, sep="")
##' grp$mean <- 0
##' grp$cov <- diag(1)
##' grp$free <- grp$param != 0
##' grp$data <- rpf.sample(500, grp=grp)
##' SitemFit(grp)
SitemFit <- function(grp, ..., method="pearson", log=TRUE, qwidth=6, qpoints=49L,
		     alt=FALSE, omit=0L, .twotier=TRUE, .parallel=TRUE) {
	if (length(list(...)) > 0) {
		stop(paste("Remaining parameters must be passed by name", deparse(list(...))))
	}

	if (missing(qwidth) && !is.null(grp$qwidth)) { qwidth <- grp$qwidth }
	if (missing(qpoints) && !is.null(grp$qpoints)) { qpoints <- grp$qpoints }

    spec <- grp$spec
    param <- grp$param
    if (ncol(param) != length(spec)) stop("Dim mismatch between param and spec")
    if (is.null(colnames(param))) stop("grp$param must have column names")

	.la <- ifelse(.parallel, mclapply, lapply)
	got <- .la(1:length(spec), function(interest) {
		free <- 0
		if (!is.null(grp$free)) free <- sum(grp$free[,interest])
		itemname <- colnames(param)[interest]
		omit1 <- omit
		if (is.list(omit)) {
			omit1 <- omit[[interest]]
		}
		ot.out <- SitemFit1(grp, itemname, free, method=method, log=log, qwidth=qwidth, qpoints=qpoints,
				    alt=alt, omit=omit1, .twotier=.twotier)
		ot.out
	})
	lapply(got, function(d1) {
		if (inherits(d1, "try-error")) stop(d1)
	})
	names(got) <- colnames(param)
	class(got) <- "summary.SitemFit"
	got
}

"+.summary.SitemFit" <- function(e1, e2) {
	e2name <- deparse(substitute(e2))
	if (!inherits(e2, "summary.SitemFit")) {
		stop("Don't know how to add ", e2name, " to a SitemFit",
		     call. = FALSE)
	}
	if (length(e1) != length(e2)) {
		stop("Cannot combine two groups with a different number of items")
	}
	if (any(names(e1) != names(e2))) {
		stop("Cannot combine two groups with a different items")
	}
	if (!all(mapply(function(i1,i2){ isTRUE(all.equal(i1$omitted, i2$omitted)) }, e1, e2))) {
		stop("Cannot combine two groups with different omitted items")
	}

	got <- mapply(function(i1, i2){
		ii <- list(orig.observed = i1$orig.observed + i2$orig.observed,
			   orig.expected = i1$orig.expected + i2$orig.expected,
			   log = i1$log,
			   method = i1$method,
			   n = i1$n + i2$n,
			   free = max(i1$free, i2$free),
			   alt = i1$alt)
		SitemFit1Internal(ii)
	}, e1, e2, SIMPLIFY=FALSE)

	class(got) <- "summary.SitemFit"
	got
}

print.summary.SitemFit <- function(x,...) {
	cat("Orlando & Thissen (2000) sum-score based item fit test\n")
	cat("  Magnitudes larger than abs(log(.01))=4.6 are significant at the p=.01 level\n\n")
	width <- max(sapply(names(x), nchar))
	fmt <- paste("%", width, "s : n = %4d, ", sep="")
	for (ix in 1:length(x)) {
		report1 <- x[[ix]]
		msg <- sprintf(fmt, names(x)[ix], report1$n)
		stat <- round(report1$statistic, 2)
		if (report1$method == "pearson") {
			msg <- paste(msg, sprintf("S-X2(%3d) = %6.2f, ", report1$df, stat), sep="")
		} else if (report1$method == "rms") {
			msg <- paste(msg, "MS=", stat, ", ", sep="")
		} else {
			stop(report1$method)
		}
		if (report1$log) {
			msg <- paste(msg, "log(p) = ", round(report1$pval, 2), sep="")
		} else {
			msg <- paste(msg, "p = ", round(report1$pval, 4), sep="")
		}
		msg <- paste(msg, "\n", sep="")
		cat(msg)
	}
}

##' Compute the ordinal gamma association statistic
##'
##' @param mat a cross tabulation matrix
##' @references
##' Agresti, A. (1990). Categorical data analysis. New York: Wiley.
##' @examples
##' # Example data from Agresti (1990, p. 21)
##' jobsat <- matrix(c(20,22,13,7,24,38,28,18,80,104,81,54,82,125,113,92), nrow=4, ncol=4)
##' ordinal.gamma(jobsat)
ordinal.gamma <- function(mat) .Call(ordinal_gamma_wrapper, mat)

# root mean squared statistic (sqrt omitted)
ms <- function(observed, expected, draws) {
  draws * sum((observed - expected)^2)
}

P.cdf.fn <- function(x, g.var, t) {
  got <- sapply(t, function (t1) {
    n <- length(g.var)
    num <- exp(1-t1) * exp(1i * t1 * sqrt(n))
    den <- pi * (t1 - 1/(1-1i*sqrt(n)))
    pterm <- prod(sqrt(1 - 2*(t1-1)*g.var/x + 2i*t1*g.var*sqrt(n)/x))
    Im(num / (den * pterm))
  })
#  print(cbind(t,got))
  got
}

##' Compute the P value that the observed and expected tables come from the same distribution
##'
##' This test is an alternative to Pearson's X^2
##' goodness-of-fit test.  In contrast to Pearson's X^2, no ad hoc cell
##' collapsing is needed to avoid an inflated false positive rate
##' in situations of sparse cell frequences.
##' The statistic rapidly converges to the Monte-Carlo estimate
##' as the number of draws increases.
##' 
##' @param observed observed matrix
##' @param expected expected matrix
##' @return The P value indicating whether the two tables come from
##' the same distribution. For example, a significant result (P <
##' alpha level) rejects the hypothesis that the two matrices are from
##' the same distribution.
##' @references Perkins, W., Tygert, M., & Ward, R. (2011). Computing
##' the confidence levels for a root-mean-square test of
##' goodness-of-fit. \emph{Applied Mathematics and Computations,
##' 217}(22), 9072-9084.
##' @examples
##' draws <- 17
##' observed <- matrix(c(.294, .176, .118, .411), nrow=2) * draws
##' expected <- matrix(c(.235, .235, .176, .353), nrow=2) * draws
##' ptw2011.gof.test(observed, expected)  # not signficiant

ptw2011.gof.test <- function(observed, expected) {
  orig.draws <- sum(observed)
  oeDiff <- abs(sum(expected) - orig.draws)
  if (is.na(oeDiff) || oeDiff > 1e-6) {
	  warning(paste("Total observed - total expected", oeDiff))
    return(NA)
  }
  if (any(c(expected)==0)) {
	  zeros <- sum(c(expected)==0)
	  warning(paste("There are", zeros, "zeros in the expected distribution.",
		     "Did you swap the observed and expected arguments"))
	  return(NA)
  }
  observed <- observed / orig.draws
  expected <- expected / orig.draws

  X <- ms(observed, expected, orig.draws)
  if (X == 0) return(1)

  n <- length(c(observed))
  D <- diag(1/c(expected))
  P <- matrix(-1/n, n,n)
  diag(P) <- 1 - 1/n
  B <- P %*% D %*% P
  g.var <- 1 / eigen(B, only.values=TRUE)$values[-n]
  
  # Eqn 8 needs n variances, but matrix B (Eqn 6) only has n-1 non-zero
  # eigenvalues. Perhaps n-1 degrees of freedom?
  
# debugging:
#  plot(function(t) P.cdf.fn(X, g.var, t), 0, 40)

  # 310 points should be good enough for 500 bins
  # Perkins, Tygert, Ward (2011, p. 10)

# If integration tolerance is too large, non-convergence can result,
# http://r.789695.n4.nabble.com/Need-help-to-understand-integrate-function-td2322093.html
  got <- try(integrate(function(t) P.cdf.fn(X, g.var, t), 0, 40, subdivisions=310L, rel.tol=1e-10), silent=TRUE)
  if (inherits(got, "try-error")) return(NA)
  p.value <- 1 - got$value
  smallest <- 6.3e-16  # approx exp(-35)
  if (p.value < smallest) p.value <- smallest
  p.value
}

CT1997Internal1 <- function(info, method) {
	observed <- info$orig.observed
	expected <- info$orig.expected
	
	s <- ordinal.gamma(observed) - ordinal.gamma(expected)
	if (!is.finite(s) || is.na(s) || s==0) s <- 1
	info <- c(info, sign=sign(s), gamma=s)

	if (method == "pearson") {
		kc <- .Call(collapse_wrapper, observed, expected)
		observed <- kc$O
		expected <- kc$E
		mask <- !is.na(expected)
		x2 <- sum((observed[mask] - expected[mask])^2 / expected[mask])
		df <- prod(dim(observed)-1) - kc$collapsed
		if (df < 1L) df <- 1L
		info <- c(info, list(statistic=x2, df=df, observed=observed, expected=expected))
	} else if (method == "lr") {
		mask <- observed > 0
		g2 <- -2 * sum(observed[mask] * log(expected[mask] / observed[mask]))
		df <- prod(dim(observed)-1)
		info <- c(info, statistic=g2, df=df)
	}
	info
}

CT1997Internal2 <- function(inames, detail) {
	cnames <- inames[-length(inames)]
	gamma <- matrix(NA, length(inames), length(cnames))
	dimnames(gamma) <- list(inames, cnames)
	raw <- matrix(NA, length(inames), length(cnames))
	dimnames(raw) <- list(inames, cnames)
	std <- matrix(NA, length(inames), length(cnames))
	dimnames(std) <- list(inames, cnames)
	pval <- matrix(NA, length(inames), length(cnames))
	dimnames(pval) <- list(inames, cnames)

	px <- 1L
	for (iter1 in 2:length(inames)) {
		for (iter2 in 1:(iter1-1)) {
			d1 <- detail[[px]]
			gamma[iter1, iter2] <- d1$gamma
			s <- d1$sign
			stat <- d1$statistic
			df <- d1$df
			raw[iter1, iter2] <- stat
			std[iter1, iter2] <- s * abs((stat - df)/sqrt(2*df))
			pval[iter1, iter2] <- s * -pchisq(stat, df, lower.tail=FALSE, log.p=TRUE)
			px <- px + 1L
		}
	}

	retobj <- list(pval=pval[-1,,drop=FALSE], std=std[-1,,drop=FALSE],
		       raw=raw[-1,,drop=FALSE], gamma=gamma[-1,,drop=FALSE], detail=detail)
	class(retobj) <- "summary.ChenThissen1997"
	retobj
}

tableWithWeights <- function(colpair, weights) {
	if (length(colpair) != 2) stop("Not a pair")
	l1 <- levels(colpair[[1]])
	l2 <- levels(colpair[[2]])
	if (1) {
		result <- .Call(fast_tableWithWeights, colpair[[1]], colpair[[2]], weights)
	} else {
		result <- matrix(0.0, length(l1), length(l2))
		if (nrow(colpair)) for (rx in 1:nrow(colpair)) {
			row <- colpair[rx,]
			ind <- sapply(row, unclass)
			w <- 1
			if (length(weights)) w <- weights[rx]
			result[ind[1], ind[2]] <- result[ind[1], ind[2]] + w
		}
	}
	dimnames(result) <- list(l1, l2)
	names(dimnames(result)) <- colnames(colpair)
	result
}

##' Computes local dependence indices for all pairs of items
##'
##' Item Factor Analysis makes two assumptions: (1) that the latent
##' distribution is reasonably approximated by the multivariate Normal
##' and (2) that items are conditionally independent. This test
##' examines the second assumption. The presence of locally dependent
##' items can inflate the precision of estimates causing a test to
##' seem more accurate than it really is.
##'
##' Statically significant entries suggest that the item pair has
##' local dependence. Since log(.01)=-4.6, an absolute magitude of 5
##' is a reasonable cut-off. Positive entries indicate that the two
##' item residuals are more correlated than expected. These items may share an
##' unaccounted for latent dimension. Consider a redesign of the items
##' or the use of testlets for scoring. Negative entries indicate that
##' the two item residuals are less correlated than expected.
##'
##' @param grp a list with the spec, param, mean, and cov describing the group
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param data data
##' @param inames a subset of items to examine
##' @param qwidth quadrature width
##' @param qpoints number of equally spaced quadrature points
##' @param method method to use to calculate P values. The default is the
##' Pearson X^2 statistic. Use "lr" for the similar likelihood ratio statistic.
##' @param .twotier whether to enable the two-tier optimization
##' @param .parallel whether to take advantage of multiple CPUs (default TRUE)
##' @return a list with raw, pval and detail. The pval matrix is a
##' lower triangular matrix of log P values with the sign
##' determined by relative association between the observed and
##' expected tables (see \code{\link{ordinal.gamma}})
##' @aliases chen.thissen.1997
##' @references Chen, W.-H. & Thissen, D. (1997). Local dependence
##' indexes for item pairs using Item Response Theory. \emph{Journal
##' of Educational and Behavioral Statistics, 22}(3), 265-289.
##'
##' Thissen, D., Steinberg, L., & Mooney, J. A. (1989). Trace lines for testlets: A use
##' of multiple-categorical-response models. \emph{Journal of Educational Measurement,
##' 26} (3), 247--260.
##'
##' Wainer, H. & Kiely, G. L. (1987). Item clusters and computerized
##' adaptive testing: A case for testlets.  \emph{Journal of
##' Educational measurement, 24}(3), 185--201.
##' @seealso \href{https://github.com/jpritikin/ifaTools}{ifaTools}
ChenThissen1997 <- function(grp, ..., data=NULL, inames=NULL, qwidth=6, qpoints=49, method="pearson",
			    .twotier=TRUE, .parallel=TRUE) {
	if (length(list(...)) > 0) {
		stop(paste("Remaining parameters must be passed by name", deparse(list(...))))
	}

  if (is.null(colnames(grp$param))) stop("Item parameter columns must be named")

  if (missing(data)) {
      data <- grp$data
  }
	if (missing(qwidth) && !is.null(grp$qwidth)) { qwidth <- grp$qwidth }
	if (missing(qpoints) && !is.null(grp$qpoints)) { qpoints <- grp$qpoints }

  if (method != "pearson" && method != "lr") stop(paste("Unknown method", method))
  if (missing(inames)) {
    inames <- colnames(grp$param)
  }
  if (length(inames) < 2) stop("At least 2 items are required")

  spec <- grp$spec
  if (length(spec) < dim(grp$param)[2]) {
      if (dim(grp$param)[2] %% length(spec) != 0) stop("Length of spec must match # of items")
      rep <- dim(grp$param)[2] %/% length(spec)
      while (rep > 1) {
          spec <- c(spec, grp$spec)
          rep <- rep - 1
      }
  }

	if (!is.data.frame(data)) {
		data <- as.data.frame(data)  #safe? TODO
	}
	dataMap <- match(colnames(grp$param), colnames(data))

	items <- match(inames, colnames(grp$param))

	# If we move this whole loop into C then we can avoid
	# repeated set up of the quadrature
	pairs <- list()
	for (iter1 in 2:length(items)) {
		for (iter2 in 1:(iter1-1)) {
			pairs[[ 1+length(pairs) ]] <- c(iter1, iter2)
		}
	}
	.la <- ifelse(.parallel, mclapply, lapply)
	detail <- .la(pairs, function(pair) {
		iter1 <- pair[1]
		iter2 <- pair[2]
		i1 <- items[iter1]
		i2 <- items[iter2]
		obpair <- data[,dataMap[c(i1,i2)]]
		rowWeight <- c()
		if (!is.null(grp[['weightColumn']])) {
			rowWeight <- data[[ grp[['weightColumn']] ]]
		}
		observed <- tableWithWeights(obpair, rowWeight)
		N <- sum(observed)

		expected <- N * pairwiseExpected(grp, c(i1, i2), qwidth, qpoints, .twotier)
		if (any(dim(observed) != dim(expected))) {
			if (dim(observed)[1] != dim(expected)[1]) {
				bad <- i1
				margin <- 1
			} else {
				bad <- i2
				margin <- 2
			}
			Eoutcomes <- dim(expected)[margin]
			lev <- dimnames(observed)[[margin]]
			stop(paste(colnames(grp$param)[bad], " has ", Eoutcomes,
				   " outcomes in the model but the data has ", length(lev), " (",
				   paste(lev, collapse=", "),")", sep=""))
		}
		dimnames(expected) <- dimnames(observed)
		info <- list(orig.observed=observed, orig.expected=expected)
		info <- CT1997Internal1(info, method)
		info
	})

	lapply(detail, function(d1) {
		if (inherits(d1, "try-error")) stop(d1)
	})

	names(detail) <- sapply(pairs, function(pair) {
		paste(inames[pair[1]], inames[pair[2]], sep=":")
	})

	retobj <- CT1997Internal2(inames, detail)
	retobj$inames <- inames
	retobj$method <- method
	retobj
}

# deprecated name
chen.thissen.1997 <- ChenThissen1997

"+.summary.ChenThissen1997" <- function(e1, e2) {
	e2name <- deparse(substitute(e2))
	if (!inherits(e2, "summary.ChenThissen1997")) {
		stop("Don't know how to add ", e2name, " to a ChenThissen1997",
		     call. = FALSE)
	}
	if (length(e1$detail) != length(e2$detail)) {
		stop("Cannot combine two groups with a different number of items")
	}
	if (any(names(e1$detail) != names(e2$detail))) {
		stop("Cannot combine two groups with a different items")
	}
	method <- e1$method
	detail <- mapply(function(i1, i2) {
		ii <- list(orig.observed = i1$orig.observed + i2$orig.observed,
			   orig.expected = i1$orig.expected + i2$orig.expected)
		ii <- CT1997Internal1(ii, method)
	}, e1$detail, e2$detail, SIMPLIFY=FALSE)

	retobj <- CT1997Internal2(e1$inames, detail)
	retobj$inames <- e1$inames
	retobj$method <- method
	retobj
}

print.summary.ChenThissen1997 <- function(x,...) {
	cat("Chen & Thissen (1997) local dependence test\n")
	cat("  Magnitudes larger than abs(log(.01))=4.6 are significant at the p=.01 level\n")
	cat("  A positive (negative) sign indicates more (less) observed correlation than expected\n\n")
	print(round(x$pval,2))
}

##' Monte-Carlo test for cross-tabulation tables
##'
##' This is for developers.
##'
##' @param ob observed table
##' @param ex expected table
##' @param trials number of Monte-Carlo trials
crosstabTest <- function(ob, ex, trials) {
	if (missing(trials)) trials <- 10000
	.Call(crosstabTest_wrapper, ob, ex, trials)
}

pairwiseExpected <- function(grp, items, qwidth=6, qpoints=49L, .twotier=FALSE) {
	.Call(pairwiseExpected_wrapper, grp, qwidth, qpoints, items - 1L, .twotier)
}

CaiHansen2012 <- function(grp, method, .twotier = FALSE) {
	.Call(CaiHansen2012_wrapper, grp, method, .twotier)
}

##' Multinomial fit test
##'
##' For degrees of freedom, we use the number of observed statistics
##' (incorrect) instead of the number of possible response patterns
##' (correct) (see Bock, Giibons, & Muraki, 1998, p. 265). This is not
##' a huge problem because this test is becomes poorly calibrated when
##' the multinomial table is sparse. For more accurate p-values, you
##' can conduct a Monte-Carlo simulation study (see examples).
##'
##' Rows with missing data are ignored.
##' 
##' The full information test is described in Bartholomew & Tzamourani
##' (1999, Section 3).
##'
##' For CFI and TLI, you must provide an independence model.
##' 
##' @param grp a list with the spec, param, mean, and cov describing the group
##' @param independenceGrp a list with the spec, param, mean, and cov describing the independence group
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param method lr (default) or pearson
##' @param log whether to report p-value in log units
##' @param .twotier whether to use the two-tier optimization (default TRUE)
##' @references Bartholomew, D. J., & Tzamourani, P. (1999). The
##' goodness-of-fit of latent trait models in attitude
##' measurement. \emph{Sociological Methods and Research, 27}(4), 525-546.
##'
##' Bock, R. D., Gibbons, R., & Muraki, E. (1988). Full-information
##' item factor analysis. \emph{Applied Psychological Measurement, 12}(3),
##' 261-280.
##' @examples
##' # Create an example IFA group
##' grp <- list(spec=list())
##' grp$spec[1:10] <- rpf.grm()
##' grp$param <- sapply(grp$spec, rpf.rparam)
##' colnames(grp$param) <- paste("i", 1:10, sep="")
##' grp$mean <- 0
##' grp$cov <- diag(1)
##' grp$uniqueFree <- sum(grp$param != 0)
##' grp$data <- rpf.sample(1000, grp=grp)
##' 
##' # Monte-Carlo simulation study
##' mcReps <- 3    # increase this to 10,000 or so
##' stat <- rep(NA, mcReps)
##' for (rx in 1:mcReps) {
##'    t1 <- grp
##'    t1$data <- rpf.sample(grp=grp)
##'    stat[rx] <- multinomialFit(t1)$statistic
##' }
##' sum(multinomialFit(grp)$statistic > stat)/mcReps   # better p-value

multinomialFit <- function(grp, independenceGrp, ..., method="lr", log=TRUE, .twotier=TRUE) {
	if (length(list(...)) > 0) {
		stop(paste("Remaining parameters must be passed by name", deparse(list(...))))
	}
	todo <- list(grp)
	if (!missing(independenceGrp)) todo <- list(grp, independenceGrp)
	todo <- lapply(todo, function(gx) {
		if (is.null(gx$weightColumn)) {
			wc <- "freq"
			gx$data <- compressDataFrame(gx$data, wc)
			gx$weightColumn <- wc
			gx$observedStats <- nrow(gx$data)
		}
		if (is.null(gx$uniqueFree)) {
			warning("Number of free parameters not available; assuming 0")
			gx$uniqueFree <- 0
		}
		if (is.null(gx$observedStats)) {
			warning("Number of observed statistics unknown; assuming the number of possible response patterns")
			gx$observedStats <- prod(sapply(gx$spec, function(sp) sp$outcomes))
		}
		gx
	})
	got <- lapply(todo, CaiHansen2012, method, .twotier)
	stat <- got[[1]]$stat
	df <- todo[[1]]$observedStats - todo[[1]]$uniqueFree
	out <- list()
	out <- c(out, list(statistic=stat, df=df))
	out$pval <- pchisq(stat, out$df, lower.tail=FALSE, log.p=log)
	out$log <- log
	out$method <- method
	out$n <- got[[1]]$n
	out$omitted <- grp$omitted
	if (length(todo) == 1) {
		fi <- computeFitStatistics(stat, df, out$n, NA, NA)
	} else {
		fi <- computeFitStatistics(stat, df, out$n,
					   got[[2]]$stat, todo[[2]]$observedStats - todo[[2]]$uniqueFree)
	}
	for (k in names(fi)) out[[k]] <- fi[[k]]
	class(out) <- "summary.multinomialFit"
	out
}

print.summary.multinomialFit <- function(x,...) {
	cat("Full information multinomial fit test\n")
	part1 <- paste("n = ", x$n, ", ", x$method, "(", x$df, ") = ", round(x$statistic, 2), sep="")
	if (x$log) {
		part2 <- paste("log(p) = ", round(x$pval,2), sep="")
	} else {
		part2 <- paste("p = ", round(x$pval,4), sep="")
	}
	cat(paste(part1, ", ", part2, sep=""), fill=TRUE)
	catFitStatistics(x)
	if (!is.null(x$omitted)) {
		cat(paste("omitted: ", paste(x$omitted, collapse=", "), "\n", sep=""))
	}
}
