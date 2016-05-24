#' Bootstrap Equating Error
#' 
#' These functions return bootstrap standard errors, bias, and RMSE of
#' equating. A summary method estimates mean and weighted mean errors over the
#' score scale. 
#' 
#' Samples are drawn of size \code{xn} and \code{yn}, with replacement, from
#' each score distribution. Form Y equivalents of each form X score are then
#' obtained using either the arguments in the equating output or those
#' provided. This process is repeated \code{reps} times. Standard errors are
#' calculated as standard deviations over replications for each score point;
#' bias is the mean equated score over replications, minus the criterion; and
#' RMSE is the square root of the squared standard error and squared bias
#' combined.
#' 
#' The bootstrap method for objects of class \dQuote{\code{equate}} is designed
#' to be called from within \code{\link{equate}}. It simply extracts the
#' necessary arguments from the equating output before bootstrapping.
#' 
#' When each element in \code{args} is a named list of equating arguments,
#' multiple equatings are performed at each replication in the bootstrapping.
#' 
#' The summary method returns a \code{data.frame} of mean standard errors,
#' bias, and rmse, and weighted and absolute means, as applicable.
#' 
#' @param x either an equating object, obtained with the \code{\link{equate}}
#' function, or a score distribution of class \dQuote{\code{\link{freqtab}}}.
#' @param xp,yp optional frequency tables replacing those equated in \code{x},
#' used for parametric bootsampling.
#' @param y score distribution of class \dQuote{\code{\link{freqtab}}}.
#' @param xn,yn integers specifying the number of scores to sample from each
#' distribution at each replication (default is the total number observed in
#' each).
#' @param reps number of bootstrap replications.
#' @param crit vector of equated scores serving as the criterion equating
#' function when calculating bootstrap bias and RMSE, both of which are
#' returned when \code{crit} is specified.
#' @param args named list of equating arguments, passed to
#' \code{\link{equate}}, specifying, e.g., the equating type and method. See
#' below for details.
#' @param eqs logical, with default \code{FALSE}, indicating whether or not the
#' matrices of equating functions (one column per replication, per equating)
#' should be returned.
#' @param object \code{bootstrap} output to be summarized.
#' @param weights vector of weights to be used in calculating weighted average
#' errors with \code{summary}, defaulting to the frequencies in
#' \code{margin(object$x)}.
#' @param subset vector indicating a subset of the score scale for which errors
#' should be summarized.
#' @param \dots further arguments passed to or from other methods.
#' @return With \code{bootstrap}, a list is returned, containing arguments
#' supplied for \code{x}, \code{y}, \code{reps}, \code{xn}, \code{yn}, and
#' \code{args}. For a single equating, the \code{mean} equating function over
#' replications and a vector of standard errors \code{se} are included,
#' along with vectors of \code{bias} and \code{rmse}, when \code{crit} is
#' provided, and a matrix of equating functions \code{eqs} when
#' \code{eqs = TRUE}. For multiple equatings, where each element of
#' \code{args} is a list of equating arguments, matrices are returned for the
#' mean functions, standard error, bias, and RMSE, and the equating functions
#' will be returned as a list of matrices. The \code{summary} method returns a
#' data frame of mean standard errors, bias, and rmse, and weighted and
#' absolute means, as applicable.
#' @author Anthony Albano \email{tony.d.albano@@gmail.com}
#' @seealso \code{\link{plot.bootstrap}}
#' @keywords methods
#' @examples
#' 
#' # Parametric bootstrapping using smoothed
#' # frequency distributions
#' set.seed(111213)
#' x <- freqtab(KBneat$x, scales = list(0:36, 0:12))
#' y <- freqtab(KBneat$y, scales = list(0:36, 0:12))
#' xp <- loglinear(x, asfreqtab = TRUE)
#' yp <- loglinear(y, asfreqtab = TRUE)
#' crit <- equate(xp, yp, "e", "c")$conc$yx
#' eqargs <- list(m.t = list(type = "m", method = "t"),
#'   l.t = list(type = "l", method = "t"))
#' bootout1 <- bootstrap(x = x, y = y, xn = 20, yn = 20,
#'   crit = crit, args = eqargs, reps = 30)
#' plot(bootout1, out = "rmse", legendplace = "top",
#'   addident = FALSE)
#' 
#' # Bootstraps for an existing equating
#' eq <- equate(x, y, type = "m", method = "t")
#' bootout2 <- bootstrap(eq, xn = 100, yn = 100,
#'   crit = crit, reps = 20)
#' summary(bootout2)
#' 
#' @export
bootstrap <- function(x, ...) UseMethod("bootstrap")

#' @describeIn bootstrap Default boostrap method for
#' \dQuote{\code{\link{freqtab}}} objects.
#' @export
bootstrap.default <- function(x, y, ...) {

	if(!is.freqtab(x) | is.freqtab(y))
		stop("'x' and 'y' must be frequency tables")
	else do.call(bootstrap.freqtab, c(list(x = x, y = y),
		list(...)))
}
		
#----------------------------------------------------------------
# Method for equate class

#' @describeIn bootstrap Method for \dQuote{\code{\link{equate}}} objects.
#' @export
bootstrap.equate <- function(x, xp = x$x, yp = x$y, ...) {
	
	dots <- list(...)
	if(is.character(xp))
		xp <- x[[xp]]
	if(is.character(yp))
		yp <- x[[yp]]
	rmnames <- c("x", "y", "yx", "concordance",
		"bootstraps", "coefficients", "synthstats",
		"xsynthetic", "ysynthetic", "xsmooth", "ysmooth",
		"points")
	args <- x[-pmatch(rmnames, names(x), nomatch = 0)]
	dots[pmatch(rmnames, names(dots), nomatch = 0)] <- NULL
	mi <- pmatch(names(dots), names(args), nomatch = 0)
	args[mi] <- dots[as.logical(mi)]
	dots <- dots[!as.logical(mi)]
	do.call(bootstrap.freqtab, c(list("x" = xp, "y" = yp),
		args, dots))
}

#----------------------------------------------------------------
# Method for freqtab class

#' @describeIn bootstrap Bootstrap method for \dQuote{\code{\link{freqtab}}}
#' objects.
#' @export
bootstrap.freqtab <- function(x, y, xn = sum(x),
	yn = sum(y), reps = 100, crit, args,
	eqs = FALSE, ...) {
	
	dots <- list(...)[names(list(...) != "")]
	if(missing(args)) {
		args <- list(dots)
		neq <- 1
		args[[1]]["verbose"] <- FALSE
	}
	else {
		neq <- length(args)
		for(i in 1:neq) {
			args[[i]][names(dots)] <- dots
			args[[i]]["verbose"] <- FALSE
		}
	}
	if(missing(y)) {
		yn <- xn
		y <- NULL
		xs <- scales(x, 1)
		ys <- scales(x, 2)
		xd <- as.data.frame(as.data.frame(x)[x > 0, 1:2])
		xp <- x[x > 0]/sum(x)
		xni <- nrow(xd)
		eqmats <- lapply(rep(NA, neq), matrix,
			nrow = length(xs), ncol = reps)
		for(i in 1:reps) {
			xi <- sample.int(xni, xn, replace = TRUE, prob = xp)
			xtemp <- freqtab(xd[xi, ], scales = list(xs, ys))
			for(j in 1:neq)
				eqmats[[j]][, i] <- do.call("equate",
					c(list(x = xtemp), args[[j]]))
		}
	}
	else {
		nx <- margins(x)
		ny <- margins(y)
		xs <- scales(x, 1:nx)
		ys <- scales(y, 1:ny)
		xd <- as.data.frame(as.data.frame(x)[x > 0, 1:nx])
		yd <- as.data.frame(as.data.frame(y)[y > 0, 1:ny])
		xp <- x[x > 0]/sum(x)
		yp <- y[y > 0]/sum(y)
		xni <- nrow(xd)
		yni <- nrow(yd)
		eqmats <- lapply(rep(NA, neq), matrix,
			nrow = length(scales(x, 1)), ncol = reps)
		for(i in 1:reps) {
			xi <- sample.int(xni, xn, replace = TRUE, prob = xp)
			xtemp <- freqtab(xd[xi, ], scales = xs)
			yi <- sample.int(yni, yn, replace = TRUE, prob = yp)
			ytemp <- freqtab(yd[yi, ], scales = ys)
			for(j in 1:neq)
				eqmats[[j]][, i] <- do.call("equate",
					c(list(x = xtemp, y = ytemp), args[[j]]))
		}
	}
	names(eqmats) <- names(args)
	out <- list(x = x, y = y, reps = reps, xn = xn, yn = yn,
		args = args, mean = sapply(eqmats, apply, 1, mean),
		se = sapply(eqmats, apply, 1, sd))
	if(!missing(crit)) {
		out$bias <- sapply(eqmats, apply, 1, mean) - crit
		out$rmse <- sqrt(out$bias^2 + out$se^2)
	}
	if(neq == 1)
		out[-(1:6)] <- lapply(out[-(1:6)], c)
	if(eqs)
		out$eqs <- if(neq == 1) eqmats[[1]] else eqmats
	out <- as.bootstrap(out)

	return(out)
}

#----------------------------------------------------------------
# Assign bootstrap class

as.bootstrap <- function(x) {
	
	class(x) <- "bootstrap"
	return(x)
}

#----------------------------------------------------------------
# Test for bootstrap class

is.bootstrap <- function(x) {
	
	return(class(x)[1] == "bootstrap")
}

#----------------------------------------------------------------
# Print method

#' @export
print.bootstrap <- function(x, ...) {
	
	nf <- length(x$args)
	cat("\nBootstrap Equating Error\n\n")
	cat("Design:", if(is.null(x$y)) "single group"
			else if(margins(x$x) == 1) "equivalent groups"
			else "nonequivalent groups", "\n\n")
	cat("Replications:", x$reps, "\n\n")
	cat("Sample Sizes: x =", paste(x$xn, "; y =", sep = ""),
		x$yn, "\n\n")
}

#----------------------------------------------------------------
# Summary method

# @describeIn bootstrap Summary method for \dQuote{\code{bootstrap}} objects.
#' @rdname bootstrap
#' @export
summary.bootstrap <- function(object, weights,
  subset, ...) {
  
  if(missing(subset))
    subset <- 1:length(scales(object$x))
  if(missing(weights))
    weights <- c(margin(object$x))[subset]/
    sum(margin(object$x)[subset])
  tempse <- cbind(object$se)[subset, , drop = FALSE]
  out <- data.frame(se = apply(tempse, 2, mean),
    w.se = apply(tempse * weights, 2, mean))
  if(!is.null(object$bias)) {
    tempbias <- cbind(object$bias)[subset, , drop = FALSE]
    out$bias <- apply(tempbias, 2, mean)
    out$a.bias <- apply(abs(tempbias), 2, mean)
    out$w.bias <- apply(tempbias * weights, 2, mean)
    out$wa.bias <- apply(abs(tempbias * weights), 2, mean)
    out$rmse <- apply(cbind(object$rmse)[subset, , drop = FALSE],
      2, mean)
    out$w.rmse <- apply(cbind(object$rmse)[subset, , drop = FALSE] *
        weights, 2, mean)
  }
  class(out) <- c("summary.bootstrap", "data.frame")
  
  return(out)
}
