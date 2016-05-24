#' Massively parallel semiparametric mixed models
#' 
#' Fits a set of semiparametric mixed models, with a common design matrix, by
#' repeated calls to \code{\link[gamm4]{gamm4}}. Only a single smooth term is
#' permitted.
#' 
#' Unlike \code{\link{semipar.mp}}, this function does not use large matrix
#' multiplications to avoid looping through model fits. Instead it performs a
#' separate call to \code{\link[gamm4]{gamm4}} to fit a semiparametric mixed
#' model for each column of \code{Y}.
#' 
#' @param Y \eqn{n \times V} response matrix.
#' @param x a vector giving the predictor upon which each column of \code{Y} is
#' regressed.
#' @param param a matrix or vector for the parametric terms in the model.
#' @param random a formula, passed to \code{\link[gamm4]{gamm4}}, specifying
#' the random effects structure in \code{\link[lme4]{lmer}} style. See the
#' example.
#' @param data.ran a required data frame containing the factors used for random
#' effects.
#' @param k number of knots.
#' @param norder order of B-splines: the default, \code{4}, gives cubic
#' B-splines.
#' @param pen.order order of the derivative penalty.
#' @param knots knot placement for the B-spline bases. The default,
#' \code{"quantile"}, gives knots at equally spaced quantiles of the data. The
#' alternative, \code{"equispaced"}, gives equally spaced knots.
#' @param store.gamm4 logical: should the \code{\link[gamm4]{gamm4}} objects to
#' be stored in the output? \code{FALSE} by default.
#' @return \item{coef}{matrix of the coefficients obtained from
#' \code{\link[gamm4]{gamm4}} looping (including both parametric and
#' nonparametric parts).} \item{bsplinecoef}{matrix of B-spline coefficients.}
#' \item{pwdf}{vector of pointwise effective degrees of freedom.}
#' \item{pwlsp}{vector of pointwise log smoothing parameters: grid values
#' maximizing the restricted likelihood at each point.} \item{B}{matrix of
#' basis function values.} \item{C}{the constraint matrix.}
#' \item{Z}{transformation matrix to impose constraints.} \item{basis}{B-spline
#' basis object, of the type created by the \pkg{fda} package; the coefficient
#' estimates are with respect to this basis.}
#' @author Yin-Hsiu Chen \email{enjoychen0701@@gmail.com} and Philip Reiss
#' \email{phil.reiss@@nyumc.org}
#' @examples
#' Y = matrix(rnorm(3000),,3)
#' x1 = rnorm(1000)
#' x2 = matrix(rnorm(2000),,2)
#' family.fac <- factor(rep(1:20,rep(50,20)))
#' person.fac <- factor(rep(rep(1:25,rep(2,25)),rep(20,50)))
#' semimix = semipar.mix.mp(Y = Y, x = x1, param = x2, random = ~ (1|a/b), 
#'           data.ran = data.frame(a = family.fac, b = person.fac))
#' @export

semipar.mix.mp <- function(Y, x, param=NULL, random, data.ran, k = 10, norder = 4, pen.order=2, knots = "quantile", store.gamm4 = FALSE) {
    if (store.gamm4) {
    	gamm4.list <- vector("list", NCOL(Y)) 
    }
    if (knots == "quantile") bs = "bq"
    else if (knots == "equispaced") bs = "be"
    coefmat <- matrix(NA, k + NCOL(param)-is.null(param), NCOL(Y))
    bsplinecoefmat <- matrix(NA, k, NCOL(Y))
    pwdf <- pwlsp <- rep(NA, NCOL(Y))
    ss = mgcv::smoothCon(mgcv::s(x, k = k, bs=bs, m=c(norder-2,pen.order)), data.frame(x = x), knots=NULL)
    B <- ss[[1]]$X 
    C <- ss[[1]]$C
    # all.equal(C, matrix(colSums(ss[[1]]$X), 1))  # TRUE
    Z <- qr.Q(qr(t(C)), complete=TRUE)[ , -1]   
    # Note that Z, as defined above, transforms the non-intercept coeffs 
    # to the original parametrization.
    basis <-  create.bspline.basis(range(x), norder = norder, nbasis = k)

    if(!is.null(param)) {
	   frame <- vector("list", length(all.vars(random)) + 3)
	   names(frame) <- c("y", "x", "param", all.vars(random))
         frame$x <- x
	   frame$param <- param
	   for(i in 1:length(all.vars(random))) frame[[3+i]] <- data.ran[, all.vars(random)[i]]
    }
    if (is.null(param)) {
         frame <- vector("list", length(all.vars(random)) + 2)
	   names(frame) <- c("y", "x", all.vars(random))
         frame$x <- x
	   for(i in 1:length(all.vars(random))) frame[[2+i]] <- data.ran[, all.vars(random)[i]]
    }

    for (j in 1:NCOL(Y)) {
    	if (!j%%40) cat("Vertex", j, "\n")
        if (is.null(param)) {
        	  frame$y <- Y[, j]
 	        gamm4.obj <-gamm4(y ~ s(x, k = k, bs=bs, m=c(norder-2,pen.order)), random=random, data=frame)
 	    }
 	    else {
   	       frame$y <- Y[, j]
             gamm4.obj <-gamm4(y ~ param + s(x, k = k, bs=bs, m=c(norder-2,pen.order)), random=random, data = frame)   
	 }
        coefmat[, j] <- gamm4.obj$gam$coefficients
        bsplinecoefmat[, j] <- coefmat[1, j] * rep(1, k) + Z %*% if (is.null(param)) coefmat[-1, j] else coefmat[-(1:(1+NCOL(param))), j]
        pwdf[j] <- sum(gamm4.obj$gam$edf)
        pwlsp[j] <- log(gamm4.obj$gam$sp)
        if (store.gamm4) gamm4.list[[j]] <- gamm4.obj
    } 
    list(coef = coefmat, bsplinecoef = bsplinecoefmat, pwdf = pwdf, pwlsp = pwlsp, B = B, C = C, 
         Z = Z, basis = basis, gamm4.list = if(store.gamm4) gamm4.list else NULL)
}
