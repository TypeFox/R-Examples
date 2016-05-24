#     This file is part of FAiR, a program to conduct Factor Analysis in R
#     Copyright 2008 Benjamin King Goodrich
#
#     FAiR is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Affero General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     FAiR is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Affero General Public License for more details.
#
#     You should have received a copy of the GNU Affero General Public License
#     along with FAiR.  If not, see <http://www.gnu.org/licenses/>.


## This file defines S4 methods for existing S3 generic functions that are hijacked

## NOTE: This file is meant to be read with 90 columns with 8 space tabs

setMethod("pairs", "FA",
function(x, ...) { # Thurstone style plot of reference structure correlations
	par(pty = "s")
	Upsilon <- loadings(x, matrix = "RS")
	if(ncol(Upsilon) < 2) {
		stop("pairs is sensible only for models with multiple factors")
	}
	rows <- nrow(Upsilon)

	bgs <- apply(!Upsilon, 1, FUN = function(x) {
			mark <- which(x)
			if(length(mark) == 0) return("white")
			else return(palette()[mark[1]])
		})
	colors <- apply(!Upsilon, 1, FUN = function(x) {
			mark <- which(x)
			if(length(mark) == 0) return("white")
			else if(length(mark) == 1) return(palette()[mark[1]])
			else return(palette()[mark[2]]) 
		})
	colors <- ifelse(colors == "white" & bgs == "white", "black", colors)

	D <- cormat(x, matrix = "PR")
	Upsilon <- rbind(Upsilon, D)
	D <- diag(D)
	Psi <- cormat(x, matrix = "RF")
	the.range <- c(min(c(Psi, Upsilon)), max(Upsilon))
	FC <- loadings(x, matrix = "FC")
	sizes <- format(colSums(FC) / sum(FC), digits = 3)
	check.fun <- function(z, current) {
			isTRUE(all.equal(z,current, check.attributes = FALSE))
		}

	UP <- function(x, y, ...) {
		abline(v = 0, col = "gray", lty = "dashed")
		abline(h = 0, col = "gray", lty = "dashed")
		xcol <- which(apply(Upsilon, 2, check.fun, current = x))
		ycol <- which(apply(Upsilon, 2, check.fun, current = y))
		abline(v = D[xcol], col = "red", lty = "dotted")
		abline(h = D[ycol], col = "red", lty = "dotted")
		points(x[1:rows], y[1:rows], col = colors, bg = bgs, pch = 21, ...)
		}

	LP <- function(x, y, ...) {
		xcol <- which(apply(Upsilon, 2, check.fun, current = x))
		ycol <- which(apply(Upsilon, 2, check.fun, current = y))
		xpoint <- Psi[xcol,ycol]
		ypoint <- sin(acos(xpoint))
		abline(h = 0, col = "gray", lty = "dashed")
		segments(x0 = 0, y0 = 0, x1 = xpoint, y1 = 0)
		segments(x0 = 0, y0 = 0, x1 = xpoint, y1 = ypoint)
		}

	diagmark <- 1
	diagmark.env <- new.env()
	environment(diagmark) <- diagmark.env
	TP <- function(x, y, labels, cex, font, envir = diagmark.env, ...) {
			diagmark <- get("diagmark", envir = envir)
			text(x, y, labels, col = diagmark, cex = 1.5 * cex)
			assign("diagmark", diagmark + 1, envir = envir)
		}

	pairs.default(Upsilon, xlim = the.range, ylim = the.range, text.panel = TP,
			labels = paste("Factor", 1:ncol(Upsilon), "\n", sizes),
			upper.panel = UP, lower.panel = LP, main = paste(
			"Lower panel: Intersection of reference axes at the origin\n",
			"Upper panel: Reference structure correlations"))
})

setMethod("fitted", "restrictions",
function(object, reduced = TRUE, standardized = TRUE) {
	Phi  <- cormat(object)
	beta <-   coef(object)
	
	chol_Phi <- chol(Phi)
	Lambda   <- chol_Phi %*% t(beta)      # weird orthogonal basis
	R <- crossprod(Lambda)                # reduced correlation matrix
	if(!reduced) diag(R) <- 1
	if(!standardized) R  <- R * tcrossprod(object@Omega@x)
	return(R)
})

setMethod("fitted", "restrictions.independent",
function(object, reduced = TRUE, standardized = TRUE) {
	R <- diag(length(object@Omega@x))
	if(!standardized) R <- R * tcrossprod(object@Omega@x)
	return(R)
})

setMethod("fitted", "restrictions.factanal",
function(object, ...) {
	stop("fitted method is not defined for 'restrections.factanal'")
})

setMethod("fitted", "restrictions.orthonormal",
function(object, reduced = TRUE, standardized = TRUE) {
	Lambda <- coef(object)
	R <- tcrossprod(Lambda)               # reduced correlation matrix
	if(!reduced) diag(R) <- 1
	if(!standardized) R  <- R * tcrossprod(object@Omega@x)
	return(R)
})

setMethod("fitted", "FA",
function(object, reduced = TRUE, standardized = TRUE) { # construct C matrix
	return(fitted(object@restrictions, reduced, standardized))
})

setMethod("residuals", "FA", 
function(object, standardized = TRUE) { # correlation or covariance residuals
	S <- model.matrix(object, standardized = standardized)
	C <- fitted(object, reduced = TRUE, standardized = standardized)
	R <- S - C
	return(R)
})

setMethod("rstandard", "FA",
function(model) { # standardized correlation residuals
	R <- residuals(model, standardized = FALSE)
	R <- R / tcrossprod(model@scale)
	return(R)
})

setMethod("weights", "FA",
function(object) { # (approximate) weights used
	if(FAiR_is.ML(object)) { ## FIXME: to Browne, MacCallum, etc. result later
		w <- 1/tcrossprod(uniquenesses(object))
	}
	else if(FAiR_is.QD(object)){
		l <- length(object@restrictions@criteria)
		middle <- formals(object@restrictions@criteria[[l]])$middle
		w <- matrix(0,  nrow = nrow(object@manifest@acov), 
				ncol = ncol(object@manifest@acov))
		w[lower.tri(w, object@manifest@diag)] <- middle
		w <- crossprod(w)
	}
	else if(object@restrictions@discrepancy == "YWLS") { # YWLS
		C <- fitted(object, reduced = TRUE, standardized = FALSE)
		communalities <- diag(C)
# 		w <- 1 - sweep(sweep(C^2, 1, communalities, FUN = "/"),
#                                   	  2, communalities, FUN = "/")
		w <- 1 - (C^2) / tcrossprod(communalities)
	}
	else stop("discrepancy function not recognized")
	return(w)
})

setMethod("influence", "FA", 
function(model) { # weights x residuals
	if(FAiR_is.QD(model)) {
		stop("influence() has not yet been worked out for this discrepancy ",
			"function")
	}
	return(residuals(model) * weights(model))
})

setMethod("df.residual", "restrictions",
function(object) { # get degrees of freedom
	if(object@model == "SEFA") {
		warning("degrees of freedom in a SEFA model is a ballpark figure until",
			" the SEFA model is actually estimated")
	}
	return(object@dof)
})

setMethod("df.residual", "FA",
function(object) { # get degrees of freedom
	return(object@restrictions@dof)
})

setMethod("deviance", "FA",
function(object) { # get value of ultimate criterion in the objective function
	if(object@restrictions@discrepancy == "YWLS") {
		warning("YWLS is not a formal discrepancy function")
	}
	fits <- FAiR_lexical_driver(object@restrictions, object@manifest, 
					sqrt(.Machine$double.eps))
	deviance <- fits[length(fits)]
	if(is.na(N <- object@manifest@n.obs)) {
		warning("deviance is not scaled by (N - 1) because N is unknown")
	}
	else deviance <- deviance * (N - 1)
	return(deviance)
})

setMethod("model.matrix", "manifest.basic",
function(object, standardized = TRUE) { # get sample covariance / correlation matrix
	if(standardized) return(object@cor)
	else             return(object@cov)
})

setMethod("model.matrix", "FA",
function(object, standardized = TRUE) { # get sample covariance / correlation matrix
	return(model.matrix(object@manifest, standardized))
})

## simulate is in GPLv3_or_later.R
