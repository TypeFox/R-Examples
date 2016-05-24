# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2    
"Krig" <- function(x, Y, cov.function = "stationary.cov", 
    lambda = NA, df = NA, GCV = FALSE, Z = NULL, cost = 1, knots = NA, 
    weights = NULL, m = 2, nstep.cv = 200, scale.type = "user", 
    x.center = rep(0, ncol(x)), x.scale = rep(1, ncol(x)), rho = NA, 
    sigma2 = NA, method = "REML", verbose = FALSE, mean.obj = NA, 
    sd.obj = NA, null.function = "Krig.null.function", wght.function = NULL, 
    offset = 0,  na.rm = TRUE, cov.args = NULL, 
    chol.args = NULL, null.args = NULL, wght.args = NULL, W = NULL, 
                   give.warnings = TRUE, ...)
# the verbose switch prints many intermediate steps as an aid in debugging.
#
{ 
    #
    # create output list
    out <- list()
    ###########################################################
    #  First series of steps simply store pieces of the passed
    #    information to the output list (i.e. the Krig object)
    ##########################################################
        out$call <- match.call()
    #   turn off warning based on options
        if( options()$warn < 0 ){
        	give.warnings<- FALSE
        }
    #
    # save covariance function as its name
    #
    if( !is.character( cov.function)){
    out$cov.function.name <- as.character(substitute(cov.function))
    }
    else{ 
    	out$cov.function.name<-cov.function
    	} 
    #
    # save null space function as its name
    #
    out$null.function.name <- as.character(substitute(null.function))
    #
    # save weight  function as its name if it is not a NULL
    #
    if (is.null(wght.function)) {
        out$wght.function.name <- NULL
    }
    else {
        out$wght.function.name <- as.character(substitute(wght.function))
    }
    out$W <- W
    if (verbose) {
        print(out$cov.function.name)
        print(out$null.function.name)
        print(out$wght.function.name)
    }
    #
    # logical to indicate if the 'C' argument is present in cov.function
    # -- a bit of esoteric R code!
    C.arg.missing <- all(names(formals(get(out$cov.function.name))) != 
        "C")
    if (C.arg.missing) 
        stop("Need to have C argument in covariance function\nsee Exp.cov.simple as an example")
    #
    # save parameters values possibly passed to the covariance function
    # also those added to call are assumed to be covariance arguments.
    if (!is.null(cov.args)) 
        out$args <- c(cov.args, list(...))
    else out$args <- list(...)
    #
    # default values for null space function
    out$null.args <- null.args
    #
    #       set degree of polynomial null space if this is default
    #       mkpoly is used so often is it helpful to include m argument
    #       by default in Krig call.
    if (out$null.function.name == "Krig.null.function") {
        out$null.args <- list(m = m)
        out$m <- m
    }
    #
    # default values for Cholesky decomposition, these are important
    # for sparse matrix decompositions used in Krig.engine.fixed.
    if (is.null(chol.args)) {
        out$chol.args <- list(pivot = FALSE)
    }
    else {
        out$chol.args <- chol.args
    }
    # additional arguments for weight matrix.
    out$wght.args <- wght.args
    #
    # the offset is the effective number of parameters used in the GCV
    # calculations -- unless this is part of an additive model this
    # is likely zero
    out$offset <- offset
    #
    # the cost is the multiplier applied to the GCV eff.df
    # 
    # lambda and df are two ways of parameterizing the smoothness
    # and are related by a monotonic function that unfortunately
    # depends on the locations of the data.
    # lambda can be used directly in the linear algebra, df
    # must be transformed to lambda numerically using the monotonic trransformation
    # sigma2 is the error variance and rho the multiplier for the covariance
    # method is how to determine lambda
    # the GCV logical forces the code to do the more elaborate decompositions
    # that faclitate estimating lambda -- even if a specific lambda value is
    # given.
    out$cost <- cost
    out$lambda <- lambda
    out$eff.df <- df
    out$sigma2 <- sigma2
    out$rho <- rho
    out$method <- method
    out$GCV <- GCV
    #
    # correlation model information
    #
    out$mean.obj <- mean.obj
    out$sd.obj <- sd.obj
    out$correlation.model <- !(is.na(mean.obj[1]) & is.na(sd.obj[1]))
    #
    # transformation info
    out$scale.type <- scale.type
    out$x.center <- x.center
    out$x.scale <- x.scale
    #
    # verbose block
    if (verbose) {
        cat("  Cov function arguments in call  ", fill = TRUE)
        print(out$args)
        cat(" covariance function used is : ", fill = TRUE)
        print(out$cov.function.name)
    }
    ###############################################################
    # Begin modifications and transformations of input information
    # note that many of these manipulations follow a strategy
    # of passing the Krig object (out) to a function and
    # then appending the information from this function to
    # the Krig object (usually also called "out").
    #In this way the Krig object  is built up
    # in steps and the process is easier to follow.
    ###############################################################
    # various checks on x and  Y including removal of NAs in Y
    # Here is an instance of adding to the Krig object
    # in this case also some onerous bookkeeping making sure arguments are consistent
    out2 <- Krig.check.xY(x, Y, Z, weights, na.rm, verbose = verbose)
    out <- c(out, out2)
    # transform to correlation model (if appropriate)
    # find replicates and collapse to means and pool variances.
    # Transform unique x locations and knots.
    if (out$correlation.model) {
        out$y <- Krig.cor.Y(out, verbose = verbose)
    }
    out2 <- Krig.transform.xY(out, knots, verbose = verbose)
    out <- c(out, out2)
    # NOTE: knots have been transformed after this step
    #############################################################
    #  Figure out what to do
    #############################################################
    #
    # this functions works through the logic of
    # what has been supplied for lambda
    out2 <- Krig.which.lambda(out)
    out[names(out2)] <- out2  
    # Make weight matrix for observations
    #    ( this is proportional to the inverse square root of obs covariance)
    #     if a weight function or W has not been passed then this is
    #     diag( out$weightsM) for W
    #     The checks represent a limitation of this model to
    #     the  WBW type decoposition and no replicate observations.
    out$nondiag.W <- (!is.null(wght.function)) | (!is.null(W))
    # Do not continue if there there is a nondiagonal weight matrix
    # and replicate observations.
    if (out$nondiag.W) {
        if (out$knot.model | out$fixed.model) {
            stop("Non diagonal weight matrix for observations
                    not supported\nwith knots or fixed lambda.")
        }
        if (!is.na(out$shat.pure.error)) {
            stop("Non diagonal weight matrix not implemented
                    with replicate locations")
        }
    }
    #  make weight matrix and its square root having passed checks
    out <- c(out, Krig.make.W(out, verbose = verbose))
    ########################################################
    #  You have reached the Engines where the actual computing happens!
    ########################################################
    #   Do the intensive linear algebra to find the solutions
    #   this is where all the heavy lifting happens.
    #
    #   Note that all the information is passed as a list
    #   including arguments to the cholesky decomposition
    #   used within Krig.engine.fixed
    #
    # The results are saved in the component matrices
    #
    # if method=='user' then just evaluate at single lambda
    #  fixed here means a fixed lambda
    #
    # For fixed lambda the decompositions with and without knots
    # are surprisingly similar and so are in one engine.
    ###########################################################
    if (out$fixed.model) {
        out$matrices <- Krig.engine.fixed(out, verbose = verbose)
    #  The trace of A matrix in fixed lambda case is not easily computed
    #  so set this to NA.
        out$eff.df <- NA
    }
    #
    # alternative are
    # matrix decompositions suitable for
    # evaluation at many lambdas to facilitate GCV/REML estimates  etc.
    #
    if (!out$fixed.model) {
        if (out$knot.model) {
            # the knot model engine
            out$matrices <- Krig.engine.knots(out, verbose = verbose)
            out$pure.ss <- out$matrices$pure.ss
        }
        else {
            # standard engine following the basic computations for thin plate splines
            out$matrices <- Krig.engine.default(out, verbose = verbose)
        }
    }
    #
    # store basic information about decompositions
    out$nt <- out$matrices$nt
    out$np <- out$matrices$np
    out$decomp <- out$matrices$decomp
    #
    # Now determine a logical vector to indicate coefficients tied to  the
    # the 'spatial drift' i.e. the fixed part of the model
    # that is not due to the Z covariates.
    # NOTE that the spatial drift coefficients must be the first columns of the
    # M matrix
    if (is.null(out$Z)) {
        out$ind.drift <- rep(TRUE, out$nt)
    }
    else {
        
        mZ <- ncol(out$ZM)
        out$ind.drift <- c(rep(TRUE, out$nt - mZ), rep(FALSE, 
            mZ))
    }
    if (verbose) {
        cat("null df: ", out$nt, "drift df: ", sum(out$ind.drift), 
            fill = TRUE)
    }
    #########################
    # End of engine block
    #########################
    #################################################
    # Do GCV and REML search over lambda if not fixed or if GCV variable is TRUE
    # gcv.Krig, not named well,  also does a search over likelihood for lambda.
    #################################################
    if (!out$fixed.model | out$GCV) {
        if (verbose) {
            cat("call to gcv.Krig", fill = TRUE)
        }
        gcv.out <- gcv.Krig(out, nstep.cv = nstep.cv, verbose = verbose, 
            cost = out$cost, offset = out$offset, give.warnings=FALSE)
        out$gcv.grid <- gcv.out$gcv.grid
        #   save a handy summary table of the search results
        out$lambda.est <- gcv.out$lambda.est
        out$warningTable<- gcv.out$warningTable
        if( verbose){
        	cat("summaries from grid search/optimization", fill=TRUE)
        	print(out$lambda.est)
        	print(out$warningTable)
        }
        if( give.warnings){
        	#NOTE: only print out grid search warning forthe method of interest.
        	printGCVWarnings( gcv.out$warningTable, method=method)
        }
          # assign the preferred lambda either from GCV/REML/MSE or the user value
        # NOTE: gcv/reml can be done but the estimate is
        # still evaluted at the passed user values of lambda (or df)
        # If df is passed need to calculate the implied lambda value
        if (out$method != "user") {
            out$lambda <- gcv.out$lambda.est[out$method, 1]
            out$eff.df <- out$lambda.est[out$method, 2]
        }
        else {
            if (!is.na(out$eff.df)) {
                out$lambda <- Krig.df.to.lambda(out$eff.df, out$matrices$D)
            }
            else {
                out$eff.df <- Krig.ftrace(out$lambda, out$matrices$D)
            }
        }
    }
    ##########################
    # end GCV/REML block
    ##########################
    #
    # Now we clean up what has happened and stuff 
    # information into output object.
    #
    ##########################################
    # find coefficients at prefered lambda
    # and evaluate the solution at observations
    ##########################################
    #   pass replicate group means -- no need to recalculate these.
    out2 <- Krig.coef(out, yM = out$yM, verbose = verbose)
    out <- c(out, out2)
    #######################################################################
    # fitted values and residuals and predicted values for full model and
    # also on the null space (fixed
    # effects). But be sure to do this at the nonmissing x's.
    ##################################################################
    out$fitted.values <- predict.Krig(out, x = out$x, Z = out$Z, 
        eval.correlation.model = FALSE)
    out$residuals <- out$y - out$fitted.values
    #
    # this is just M%*%d  note use of do.call using function name
    Tmatrix <- do.call(out$null.function.name, c(out$null.args, 
        list(x = out$x, Z = out$Z)))
    out$fitted.values.null <- as.matrix(Tmatrix) %*% out$d
    #
    # verbose block
    if (verbose) {
        cat("residuals", out$residuals, fill = TRUE)
    }
    #
    # find various estimates of sigma and rho
    out2 <- Krig.parameters(out)
    out <- c(out, out2)
    ################################################
    # assign the 'best' model as a default choice
    # either use the user supplied values or the results from
    # optimization
    ################################################
    passed.sigma2 <- (!is.na(out$sigma2))
    if (out$method == "user" & passed.sigma2) {
        out$best.model <- c(out$lambda, out$sigma2, out$rho)
    }
    else {
        # in this case lambda is from opt. or supplied by user
        out$best.model <- c(out$lambda, out$shat.MLE^2, out$rhohat)
    }
    # Note: values in best.model are used in subsquent functions as the choice
    # for these parameters!
    # set class
    class(out) <- c("Krig")
    return(out)
}


# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

Krig.check.xY <- function(x, Y, Z, weights, na.rm, 
    verbose = FALSE) {
    #
    # check for missing values in Y or X.
    #
    # save logical indicating where there are NA's
    # and check for NA's
    #
    ind <- is.na(Y)
    if (any(ind) & !na.rm) {
        stop("Need to remove missing values or use: na.rm=TRUE in the call")
    }
    #
    # coerce x to be a matrix
    x <- as.matrix(x)
    #
    # coerce Y to be a vector
    #
    Y <- as.matrix(Y)
    if (ncol(Y) != 1) {
        stop("Krig can handle matrix Y data")
    }
    #
    #default weights ( reciprocal variance of errors).
    #
    if (is.null(weights)) 
        weights <- rep(1, length(Y))
    #
    # check that dimensions agree
    #
    if (length(Y) != nrow(x)) {
        stop(" length of y and number of rows of x differ")
    }
    if (length(Y) != length(weights)) {
        stop(" length of y and weights differ")
    }
    #  if Z is not NULL coerce to be  a matrix
    # and check  # of rows
    if (verbose) {
        print(Z)
    }
    if (!is.null(Z)) {
        if (!is.matrix(Z)) {
            Z <- matrix(c(Z), ncol = 1)
        }
        if (length(Y) != nrow(Z)) {
            stop(" length of y and number of rows of Z differ")
        }
    }
    # if NAs can be removed then remove them and warn the user
    if (na.rm) {
        ind <- is.na(Y)
        if(all(ind)){
        	stop("Oops! All Y values are missing!")
        }
        if (any(ind)) {
            Y <- Y[!ind]
            x <- as.matrix(x[!ind, ])
            if (!is.null(Z)) {
                Z <- Z[!ind, ]
            }
            weights <- weights[!ind]
        }
    }
    #
    # check for NA's in x matrix -- there should not be any !
    if (any(c(is.na(x)))) {
        stop(" NA's in x matrix")
    }
    #
    # check for NA's in Z matrix
    if (!is.null(Z)) {
        if (any(c(is.na(Z)))) {
            stop(" NA's in Z matrix")
        }
    }
    #
    # verbose block
    if (verbose) {
        cat("Y:", fill = TRUE)
        print(Y)
        cat("x:", fill = TRUE)
        print(x)
        cat("weights:", fill = TRUE)
        cat(weights, fill = TRUE)
    }
    #
    # save x, weights  and Y w/o NAs
    N <- length(Y)
    return(list(N = N, y = Y, x = x, weights = weights, Z = Z, 
        NA.ind = ind))
}

"Krig.coef" <- function(out, lambda = out$lambda, 
    y = NULL, yM = NULL, verbose = FALSE) {
    #
    # NOTE default value of lambda used from Krig object.
    #
    # Determine whether to collapse onto means of replicates ( using y)
    # if the data has been passed use as the replicate means (yM) use that.
    # If both y and YM are null then just use out$yM
    # For readability of this function, all this tortured logic happens in
    #  Krig.ynew.
    #
    out2 <- Krig.ynew(out, y, yM)
    temp.yM <- out2$yM
    nt <- out$nt
    np <- out$np
    ndata <- ncol(temp.yM)
    u <- NA
    call.name <- out$cov.function.name
    if (verbose) {
        cat("dimension of yM in Krig.coef", fill = TRUE)
        print(dim(temp.yM))
    }
    #
    #   case when knots= unqiue x's
    # any lambda
    #
    if (out$decomp == "WBW") {
        # pad u with zeroes that corresond to null space basis functions
        # this makes it compatible with the DR decomposition.
        u <- rbind(matrix(0, nrow = out$nt, ncol = ndata), t(out$matrices$V) %*% 
            qr.q2ty(out$matrices$qr.T, out$W2 %d*% temp.yM))
        #
        #old code   beta <- out$matrices$G %*% ((1/(1 + lambda * out$matrices$D))%d*%u)
        #
        ind <- (nt + 1):np
        D2 <- out$matrices$D[ind]
        #
        # note use of efficient diagonal multiply in next line
        temp2 <- (D2/(1 + lambda * D2)) %d*% u[ind, ]
        beta2 <- out$matrices$V %*% temp2
        temp.c <- rbind(matrix(0, nrow = nt, ncol = ndata), beta2)
        temp.c <- qr.qy(out$matrices$qr.T, temp.c)
        temp.c <- out$W2 %d*% temp.c
        temp <- temp.yM - do.call(call.name, c(out$args, list(x1 = out$knots, 
            x2 = out$knots, C = temp.c)))
        temp <- out$W2 %d*% temp
        temp.d <- qr.coef(out$matrices$qr.T, temp)
    }
    #
    # case with knots
    # any lambda
    #
    if (out$decomp == "DR") {
        # X is the monster matrix ...  X = [ M | K]
        X <- cbind(do.call(out$null.function.name, c(out$null.args, 
            list(x = out$xM, Z = out$ZM))), do.call(call.name, 
            c(out$args, list(x1 = out$xM, x2 = out$knots))))
        u <- t(out$matrices$G) %*% t(X) %*% (out$weightsM %d*% 
            temp.yM)
        beta <- out$matrices$G %*% ((1/(1 + lambda * out$matrices$D)) %d*% 
            u)
        temp.d <- beta[1:nt, ]
        temp.c <- beta[(nt + 1):np, ]
        temp <- X %*% out$matrices$G %*% u
        temp <- sum(out$weightsM * (temp.yM - temp)^2)
        #### ????
        out2$pure.ss <- temp + out2$pure.ss
    }
    #
    # fixed lambda knots == unique x's
    #
    if (out$decomp == "cholesky") {
        if (lambda != out$matrices$lambda) {
            stop("New lambda can not be used with cholesky decomposition")
        }
        Tmatrix <- do.call(out$null.function.name, c(out$null.args, 
            list(x = out$knots, Z = out$ZM)))
        temp.d <- qr.coef(out$matrices$qr.VT, forwardsolve(out$matrices$Mc, 
            transpose = TRUE, temp.yM, upper.tri = TRUE))
        temp.c <- forwardsolve(out$matrices$Mc, transpose = TRUE, 
            temp.yM - Tmatrix %*% temp.d, upper.tri = TRUE)
        temp.c <- backsolve(out$matrices$Mc, temp.c)
    }
    #
    # fixed lambda with knots
    #
    if (out$decomp == "cholesky.knots") {
        if (lambda != out$matrices$lambda) {
            stop("New lambda can not be used with cholesky decomposition")
        }
        # form K matrix
        K <- do.call(call.name, c(out$args, list(x1 = out$xM, 
            x2 = out$knots)))
        Tmatrix <- do.call(out$null.function.name, c(out$null.args, 
            list(x = out$xM, Z = out$ZM)))
        wY <- out$weightsM * temp.yM
        temp0 <- t(K) %*% (out$weightsM * Tmatrix)
        temp1 <- forwardsolve(out$matrices$Mc, temp0, transpose = TRUE, 
            upper.tri = TRUE)
        qr.Treg <- qr(t(Tmatrix) %*% (out$weightsM * Tmatrix) - 
            t(temp1) %*% temp1)
        temp0 <- t(K) %*% wY
        temp3 <- t(Tmatrix) %*% wY - t(temp1) %*% forwardsolve(out$matrices$Mc, 
            temp0, transpose = TRUE, upper.tri = TRUE)
        temp.d <- qr.coef(qr.Treg, temp3)
        temp1 <- t(K) %*% (wY - out$weightsM * (Tmatrix) %*% 
            temp.d)
        temp.c <- forwardsolve(out$matrices$Mc, transpose = TRUE, 
            temp1, upper.tri = TRUE)
        temp.c <- backsolve(out$matrices$Mc, temp.c)
    }
    return(list(c = temp.c, d = temp.d, shat.rep = out2$shat.rep, 
        shat.pure.error = out2$shat.pure.error, pure.ss = out2$pure.ss))
}

Krig.cor.Y <- function(obj, verbose = FALSE) {
    # subtract mean
    if (!is.na(obj$mean.obj[1])) {
        Y <- obj$y - predict(obj$mean.obj, obj$x)
    }
    # divide by sd
    if (!is.na(obj$sd.obj[1])) {
        Y <- Y/predict(obj$sd.obj, obj$x)
    }
    Y
}

Krig.Amatrix <- function(object, x0 = object$x, lambda = NULL, 
    eval.correlation.model = FALSE, ...) {
    if (is.null(lambda)) {
        lambda <- object$lambda
    }
    M <- nrow(object$xM)
    N <- nrow(x0)
    # create output matrix
    out <- matrix(NA, N, M)
    #
    # loop through unique data locations predicting response
    # using unit vector
    # NOTE that the y vector has already been collapsed onto means.
    #
    for (k in 1:M) {
        ytemp <- rep(0, M)
        ytemp[k] <- 1
        out[, k] <- predict(object, x = x0, yM = ytemp, lambda = lambda, 
            eval.correlation.model = eval.correlation.model, 
            ...)
    }
    return(out)
}
"Krig.df.to.lambda" <- function(df, D, guess = 1, 
    tol = 1e-05) {
    if (is.list(D)) {
        D <- D$matrices$D
    }
    if (is.na(df)) 
        return(NA)
    if (df < sum(D == 0)) {
        warning("df too small to match with a lambda value")
        return(NA)
    }
    if (df > length(D)) {
        warning(" df too large to match a lambda value")
        return(NA)
    }
    l1 <- guess
    for (k in 1:25) {
        tr <- sum(1/(1 + l1 * D))
        if (tr <= df) 
            break
        l1 <- l1 * 4
    }
    l2 <- guess
    for (k in 1:25) {
        tr <- sum(1/(1 + l2 * D))
        if (tr >= df) 
            break
        l2 <- l2/4
    }
    info <- list(D = D, df = df, N = length(D))
    out <- bisection.search(log(l1), log(l2), Krig.fdf, tol = tol, 
        f.extra = info)$x
    +exp(out)
}

"Krig.engine.default" <- function(out, verbose = FALSE) {
    #
    # matrix decompositions for computing estimate
    #
    # Computational outline:( '.' is used for subscript)
    #
    # The form of the estimate is
    #    fhat(x) = sum phi.j(x) d.j  + sum psi.k(x) c.k
    #
    # the {phi.j} are the fixed part of the model usually low order polynomials
    # and is also referred to as spatial drift.
    #
    # the {psi.k} are the covariance functions evaluated at the unique observation
    # locations or 'knots'.  If xM.k is the kth unique location psi.k(x)= k(x, xM.k)
    # xM is also out$knots in the code below.
    #
    # the goal is find decompositions that facilitate rapid solution for
    # the vectors d and c. The eigen approach below was identified by
    # Wahba, Bates Wendelberger and is stable even for near colinear covariance
    # matrices.
    # This function does the main computations leading to the matrix decompositions.
    # With these decompositions the coefficients of the solution are found in
    # Krig.coef and the GCV and REML functions in Krig.gcv.
    #
    #  First is an outline calculations with equal weights
    #  T the fixed effects regression matrix  T.ij = phi.j(xM.i)
    #  K the covariance matrix for the unique locations
    # From the spline literature the solution solves the well known system
    # of two eqautions:
    #    -K( yM - Td - Kc) + lambda *Kc = 0
    #                 -T^t ( yM-Td -Kc) = 0
    #
    # Mulitple through by K inverse and substitute, these are equivalent to
    #
    #  -1-   -( yM- Td - Kc) + lambda c = 0
    #  -2-                        T^t c = 0
    #
    #
    #  A QR decomposition is done for   T= (Q.1,Q.2)R
    #   by definition  Q.2^T T =0
    #
    #  equation  -2- can be thought of as a constraint
    # with  c= Q.2 beta2
    # substitute in  -1-  and multiply through by Q.2^T
    #
    #      -Q.2^T yM  + Q.2^T K Q.2 beta2  + lambda beta2 = 0
    #
    #   Solving
    #   beta2 = {Q.2^T K Q.2 + lambda I )^ {-1} Q.2^T yM
    #
    # and so one sloves this linear system for beta2 and then uses
    #     c= Q.2 beta2
    #   to determine c.
    #
    #  eigenvalues and eigenvectors are found for M= Q.2^T K Q.2
    #     M = V diag(eta) V^T
    #  and these facilitate solving this system efficiently for
    #  many different values of lambda.
    #  create eigenvectors, D = (0, 1/eta)
    #  and G= ( 0,0) %*% diag(D)
    #         ( 0,V)
    # so that
    #
    #          beta2 = G%*% ( 1/( 1+ lambda D)) %*% u
    # with
    #
    #          u = (0, V Q.2^T W2 yM)
    #
    # Throughout keep in mind that M has smaller dimension than G due to
    # handling the null space.
    #
    # Now solve for d.
    #
    # From -1-  Td = yM - Kc - lambda c
    #      (Q.1^T) Td =  (Q.1^T) ( yM- Kc)
    #
    #   ( lambda c is zero by -2-)
    #
    #   so Rd = (Q.1^T) ( yM- Kc)
    # use qr functions to solve triangular system in R to find d.
    #
    #----------------------------------------------------------------------
    # What about errors with a general precision matrix, W?
    #
    # This is an important case because with replicated observations the
    # problem will simplify into a smoothing problem with the replicate group
    # means and unequal measurement error variances.
    #
    # the equations to solve are
    #     -KW( yM - Td - Kc) + lambda *Kc = 0
    #     -T^t W( yM-Td -Kc) =0
    #
    # Multiple through by K inverse and substitute, these are equivalent to
    #
    #  -1b-      -W( yM- Td - Kc) + lambda c = 0
    #  -2b-      (WT)^t c = 0
    #
    # Let W2 be the symmetric square root of W,  W= W2%*% W2
    # and W2.i be the inverse of W2.
    #
    #  -1c-      -(  W2 yM - W2 T d - (W2 K W2) W2.ic) + lambda W2.i c = 0
    #  -2c-      (W2T)^t  W2c = 0
    Tmatrix <- do.call(out$null.function.name, c(out$null.args, 
        list(x = out$xM, Z = out$ZM)))
    if (verbose) {
        cat(" Model Matrix: spatial drift and Z", fill = TRUE)
        print(Tmatrix)
    }
    # Tmatrix premultiplied by sqrt of wieghts
    Tmatrix <- out$W2 %d*% Tmatrix
    qr.T <- qr(Tmatrix)
    if( qr.T$rank < ncol( Tmatrix)){
      stop("Regression matrix for fixed part of model is colinear")}
    #
    #verbose block
    if (verbose) {
        cat("first 5 rows of qr.T$qr", fill = TRUE)
        print(qr.T$qr[1:5, ])
    }
    #
    # find  Q_2 K Q_2^T  where K is the covariance matrix at the knot points
    #
    tempM <- t(out$W2 %d*% do.call(out$cov.function.name, c(out$args, 
        list(x1 = out$knots, x2 = out$knots))))
    tempM <- out$W2 %d*% tempM
    tempM <- qr.yq2(qr.T, tempM)
    tempM <- qr.q2ty(qr.T, tempM)
    np <- nrow(out$knots)
    nt <- (qr.T$rank)
    if (verbose) {
        cat("np, nt", np, nt, fill = TRUE)
    }
    #
    # Full set of decompositions for
    # estimator for nonzero lambda
    tempM <- eigen(tempM, symmetric = TRUE)
    D <- c(rep(0, nt), 1/tempM$values)
    #
    # verbose block
    if (verbose) {
        cat("eigen values:", fill = TRUE)
        print(D)
    }
    #
    # Find the transformed data vector used to
    # evaluate the solution, GCV, REML  at different lambdas
    #
    
    u <- c(rep(0, nt), t(tempM$vectors) %*% qr.q2ty(qr.T, c(out$W2 %d*% 
        out$yM)))
    if (verbose) {
        cat("u vector:", fill = TRUE)
        print(u)
    }
    #
    #
    return(list(D = D, qr.T = qr.T, decomp = "WBW", V = tempM$vectors, 
        u = u, nt = nt, np = np))
}

"Krig.engine.fixed" <- function(out, verbose = FALSE, 
    lambda = NA) {
    #
    # Model:
    #     Y_k=  f_k + e_k
    #  var( e_k) = sigma^2/W_k
    #
    #   f= Td + h
    #    T is often a low order polynomial
    #   E(h)=0    cov( h)= rho *K
    #
    # let M = (lambda W^{-1} + K)
    # the surface estimate  depends on coefficient vectors d and c
    #    The implementation in Krig/fields is that K are the
    #    cross covariances among the observation locations and the knot locations
    #    H is the covariance among the knot locations.
    #    Thus if knot locs == obs locs we have the obvious collapse to
    #    the simpler form for M above.
    #
    #   With M in hand ...
    #
    #   set
    #   d =  [(T)^t M^{-1} (T)]^{-1} (T)^t M^{-1} Y
    #  this is just the generalized LS estimate for d
    #
    #   lambda= sigma**2/rho
    #  the estimate for c is
    #   c=  M^{-1}(y - Td)
    #
    # This particular numerical strategy takes advantage of
    # fast Cholesky factorizations for positive definite matrices
    # and also provides a seamless framework for sparse matrix implementations
    #
    if (is.na(lambda)) 
        lambda <- out$lambda
    call.name <- out$cov.function.name
    if (!out$knot.model) {
        ####################################################
        # case of knot locs == obs locs  out$knots == out$xM
        ####################################################
        # create T matrix
        Tmatrix <- do.call(out$null.function.name, c(out$null.args, 
            list(x = out$knots, Z = out$ZM)))
        if (verbose) {
            cat("Tmatrix:", fill = TRUE)
            print(Tmatrix)
        }
        np <- nrow(out$knots)
        nt <- ncol(Tmatrix)
        # form K
        tempM <- do.call(call.name, c(out$args, list(x1 = out$knots, 
            x2 = out$knots)))
        # form M
        diag(tempM) <- (lambda/out$weightsM) + diag(tempM)
        #
        # find cholesky factor
        #  tempM = t(Mc)%*% Mc
        #  V=  Mc^{-T}
        # call cholesky but also add in the args supplied in Krig object.
        Mc <- do.call("chol", c(list(x = tempM), out$chol.args))
        VT <- forwardsolve(Mc, x = Tmatrix, transpose = TRUE, 
            upper.tri = TRUE)
        qr.VT <- qr(VT)
        # find GLS covariance matrix of null space parameters.
        Rinv <- solve(qr.R(qr.VT))
        Omega <- Rinv %*% t(Rinv)
        #
        # now do generalized least squares for d
        # and then find c.
        d.coef <- qr.coef(qr.VT, forwardsolve(Mc, transpose = TRUE, 
            out$yM, upper.tri = TRUE))
        if (verbose) {
            print(d.coef)
        }
        c.coef <- forwardsolve(Mc, transpose = TRUE, out$yM - 
            Tmatrix %*% d.coef, upper.tri = TRUE)
        c.coef <- backsolve(Mc, c.coef)
        # return all the goodies,  include lambda as a check because
        # results are meaningless for other values of lambda
        return(list(qr.VT = qr.VT, d = c(d.coef), c = c(c.coef), 
            Mc = Mc, decomp = "cholesky", nt = nt, np = np, lambda.fixed = lambda, 
            Omega = Omega))
    }
    else {
        ####################################################
        # case of knot locs != obs locs
        ####################################################
        # create weighted T matrix
        Tmatrix <- do.call(out$null.function.name, c(out$null.args, 
            list(x = out$xM, Z = out$ZM)))
        nt <- ncol(Tmatrix)
        np <- nrow(out$knots) + nt
        # form H
        H <- do.call(call.name, c(out$args, list(x1 = out$knots, 
            x2 = out$knots)))
        # form K matrix
        K <- do.call(call.name, c(out$args, list(x1 = out$xM, 
            x2 = out$knots)))
        #
        Mc <- do.call("chol", c(list(x = t(K) %*% (out$weightsM * 
            K) + lambda * H), out$chol.args))
        # weighted Y
        wY <- out$weightsM * out$yM
        temp0 <- t(K) %*% (out$weightsM * Tmatrix)
        temp1 <- forwardsolve(Mc, temp0, transpose = TRUE, upper.tri = TRUE)
        qr.Treg <- qr(t(Tmatrix) %*% (out$weightsM * Tmatrix) - 
            t(temp1) %*% temp1)
        temp0 <- t(K) %*% wY
        temp3 <- t(Tmatrix) %*% wY - t(temp1) %*% forwardsolve(Mc, 
            temp0, transpose = TRUE, upper.tri = TRUE)
        d.coef <- qr.coef(qr.Treg, temp3)
        temp1 <- t(K) %*% (wY - out$weightsM * (Tmatrix) %*% 
            d.coef)
        c.coef <- forwardsolve(Mc, transpose = TRUE, temp1, upper.tri = TRUE)
        c.coef <- backsolve(Mc, c.coef)
        list(qr.Treg = qr.Treg, d = c(d.coef), c = c(c.coef), 
            Mc = Mc, decomp = "cholesky.knots", nt = nt, np = np, 
            lambda.fixed = lambda, Omega = NA)
    }
    #
    # should not get here.
    #
}

"Krig.engine.knots" <- function(out, verbose = FALSE) {
    #
    # matrix decompostions for computing estimate when
    # knots are present
    # QR decomposition of null space regression matrix
    Tmatrix <- do.call(out$null.function.name, c(out$null.args, 
        list(x = out$xM, Z = out$ZM)))
    qr.T <- qr(c(sqrt(out$weightsM)) * Tmatrix)
    nt <- ncol(Tmatrix)
    np <- nrow(out$knots) + nt
    if (verbose) {
        cat(nt, np, fill = TRUE)
    }
    # H is the penalty matrix in the ridge regression format
    # first part is zero because no penalty on part of estimator
    # spanned by T matrix
    H <- matrix(0, ncol = np, nrow = np)
    H[(nt + 1):np, (nt + 1):np] <- do.call(out$cov.function.name, 
        c(out$args, list(x1 = out$knots, x2 = out$knots)))
    # X is the monster ...
    X <- cbind(do.call(out$null.function.name, c(out$null.args, 
        list(x = out$xM, Z = out$ZM))), do.call(out$cov.function.name, 
        c(out$args, list(x1 = out$xM, x2 = out$knots))))
    if (verbose) {
        cat("first lines of X", fill = TRUE)
        print(X[1:5, ])
    }
    #    sqrt(weightsM) * X
    XTwX <- t(X * out$weightsM) %*% X
    #
    # then  B= G(I-D)G^T
    # New version of diagonalize may be more stable
    out2 <- fields.diagonalize2((XTwX), H)
    D <- out2$D
    if (verbose) {
        cat("D;", fill = TRUE)
        cat(out2$D, fill = TRUE)
    }
    #
    #  G should satisfy:
    #     t(G) %*% XTwX %*%G = I and  t(G)%*%H%*%G = D
    #
    #     and
    #      solve( XtwX + lambda H) =  G%*%diag( 1/(1+ lambda*D))%*%t(G)
    #
    
    #  save XG to avoid an extra multiplication.
    XG <- X %*% out2$G
    
    u <- t(XG) %*% (out$weightsM * out$yM)
    #
    # adjust pure sum of squares to be that due to replicates
    # plus that due to fitting all the basis functions without
    # any smoothing. This will be the part of the RSS that does not
    # change as lambda is varied ( see e.g. gcv.Krig)
    #
    pure.ss <- sum(out$weightsM * (out$yM - XG %*% u)^2) + out$pure.ss
    if (verbose) {
        cat("total pure.ss from reps, reps + knots ", fill = TRUE)
        print(out$pure.ss)
        print(pure.ss)
    }
    
    #
    # in this form  the solution is (d,c)= G( I + lambda D)^-1 u
    # fitted.values = X ( d,c)
    #
    # output list
    # last D eigenvalues are zero due to null space of penalty
    # OLD code:    D[(np - nt + 1):np] <- 0
    # this should be enforced to machine precision from diagonalization.
    
    
    list(u = u, D = D, G = out2$G, qr.T = qr.T, decomp = "DR", 
        nt = nt, np = np, pure.ss = pure.ss)
}

"Krig.fdf" <- function(llam, info) {
    sum(1/(1 + exp(llam) * info$D)) - info$df
}

"Krig.fgcv" <- function(lam, obj) {
    #
    # GCV that is leave-one-group out
    #
    lD <- obj$matrices$D * lam
    RSS <- sum(((obj$matrices$u * lD)/(1 + lD))^2)
    MSE <- RSS/length(lD)
    if ((obj$N - length(lD)) > 0) {
        MSE <- MSE + obj$pure.ss/(obj$N - length(lD))
    }
    trA <- sum(1/(1 + lD))
    den <- (1 - (obj$cost * (trA - obj$nt - obj$offset) + obj$nt)/length(lD))
    # If the denominator is negative then flag this as a bogus case
    # by making the GCV function 'infinity'
    #
    ifelse(den > 0, MSE/den^2, 1e20)
}

"Krig.fgcv.model" <- function(lam, obj) {
    lD <- obj$matrices$D * lam
    MSE <- sum(((obj$matrices$u * lD)/(1 + lD))^2)/length(lD)
    trA <- sum(1/(1 + lD))
    den <- (1 - (obj$cost * (trA - obj$nt - obj$offset) + obj$nt)/length(lD))
    ifelse(den > 0, obj$shat.pure.error^2 + MSE/den^2, 1e20)
}

"Krig.fgcv.one" <- function(lam, obj) {
    lD <- obj$matrices$D * lam
    RSS <- obj$pure.ss + sum(((obj$matrices$u * lD)/(1 + lD))^2)
    trA <- sum(1/(1 + lD))
    den <- 1 - (obj$cost * (trA - obj$nt - obj$offset) + obj$nt)/obj$N
    # If the denominator is negative then flag this as a bogus case
    # by making the GCV function 'infinity'
    #
    ifelse(den > 0, (RSS/obj$N)/den^2, 1e+20)
}

"Krig.flplike" <- function(lambda, obj) {
    #  - log profile likelihood for lambda
    # See section 3.4 from Nychka  Spatial Processes as Smoothers paper.
    # for equation and derivation
    D2 <- obj$matrices$D[obj$matrices$D > 0]
    u2 <- obj$matrices$u[obj$matrices$D > 0]
    lD <- D2 * lambda
    N2 <- length(D2)
    # MLE estimate of rho for fixed lambda
    rho.MLE <- (sum((D2 * (u2)^2)/(1 + lD)))/N2
    #
    # ln determinant of    K + lambda*WI
    lnDetCov <- -sum(log(D2/(1 + lD)))
    
    -1 * (-N2/2 - log(2 * pi) * (N2/2) - (N2/2) * log(rho.MLE) - 
        (1/2) * lnDetCov)
      
    
}

"Krig.fs2hat" <- function(lam, obj) {
    lD  <- obj$matrices$D * lam
    RSS <- obj$pure.ss + sum(((obj$matrices$u * lD)/(1 + lD))^2)
    den <- obj$N - (sum(1/(1 + lD)) + obj$offset)
    if (den < 0) {
        return(NA)
    }
    else {
        RSS/(den)
    }
}

"Krig.ftrace" <- function(lam, D) {
    sum(1/(1 + lam * D))
}

"Krig.make.W" <- function(out, verbose = FALSE) {
    if (verbose) {
        cat("W", fill = TRUE)
        print(out$W)
    }
    if (out$nondiag.W) {
        #
        # create W from scratch or grab it from passed object
        if (is.null(out$W)) {
            if (verbose) {
                print(out$wght.function.name)
            }
            W <- do.call(out$wght.function.name, c(list(x = out$xM), 
                out$wght.args))
            #       adjust W based on diagonal weight terms
            #
            W <- sqrt(out$weightsM) * t(sqrt(out$weightsM) * 
                W)
        }
        else {
            W <- out$W
        }
        #
        # symmetric square root
        temp <- eigen(W, symmetric = TRUE)
        W2 <- temp$vectors %*% diag(sqrt(temp$values)) %*% t(temp$vectors)
        return(list(W = W, W2 = W2))
    }
    else {
        #
        #  These are created only for use with default method to stay
        #   consistent with nondiagonal elements.
        if (out$fixed.model) {
            return(list(W = NULL, W2 = NULL))
        }
        else {
            return(list(W = out$weightsM, W2 = sqrt(out$weightsM)))
        }
    }
}

"Krig.make.Wi" <- function(out, verbose = FALSE) {
    #
    # If a weight matrix has been passed use it.
    #
    # Note that in either case the weight matrix assumes that
    # replicate observations have been collapses to the means.
    #
    if (out$nondiag.W) {
        temp <- eigen(out$W, symmetric = TRUE)
        Wi <- temp$vectors %*% diag(1/(temp$values)) %*% t(temp$vectors)
        W2i <- temp$vectors %*% diag(1/sqrt(temp$values)) %*% 
            t(temp$vectors)
        return(list(Wi = Wi, W2i = W2i))
    }
    else {
        #
        #  These are created only for use with default method to stay
        # consistent with nondiagonal elements.
        return(list(Wi = 1/out$weightsM, W2i = 1/sqrt(out$weightsM)))
    }
}

"Krig.make.u" <- function(out, y = NULL, yM = NULL, 
    verbose = FALSE) {
    #
    # Determine whether to collapse onto means of replicates ( using y)
    # if the data has been passed use as the replicate means (yM) use that.
    # If both y and YM are null then just use out$yM
    # For readability of this function, all this tortured logic happens in
    #  Krig.ynew.
    #
    out2 <- Krig.ynew(out, y, yM)
    temp.yM <- out2$yM
    nt <- out$nt
    np <- out$np
    ndata <- ncol(temp.yM)
    u <- NA
    call.name <- out$cov.function.name
    if (verbose) {
        cat("dimension of yM in Krig.coef", fill = TRUE)
        print(dim(temp.yM))
    }
    #
    #   case when knots= unqiue x's
    # any lambda
    #
    if (out$decomp == "WBW") {
        # pad u with zeroes that corresond to null space basis functions
        # this makes it compatible with the DR decomposition.
        u <- rbind(matrix(0, nrow = out$nt, ncol = ndata), t(out$matrices$V) %*% 
            qr.q2ty(out$matrices$qr.T, out$W2 %d*% temp.yM))
    }
    #
    # case with knots
    # any lambda
    #
    if (out$decomp == "DR") {
        # X is the monster matrix ...  X = [ M | K]
        X <- cbind(do.call(out$null.function.name, c(out$null.args, 
            list(x = out$xM, Z = out$ZM))), do.call(call.name, 
            c(out$args, list(x1 = out$xM, x2 = out$knots))))
        u <- t(out$matrices$G) %*% t(X) %*% (out$weightsM %d*% 
            temp.yM)
    }
    return(list(u = u, shat.rep = out2$shat.rep, shat.pure.error = out2$shat.pure.error, 
        pure.ss = out2$pure.ss))
}

Krig.null.function <- function(x, Z = NULL, drop.Z = FALSE, 
    m) {
    # default function to create matrix for fixed part of model
    #  x, Z, and drop.Z are required
    #  Note that the degree of the polynomial is by convention (m-1)
    # returned matrix must have the columns from Z last!
    #
    if (is.null(Z) | drop.Z) {
        return(fields.mkpoly(x, m = m))
    }
    else {
        return(cbind(fields.mkpoly(x, m = m), Z))
    }
}

"Krig.parameters" <- function(obj, mle.calc = obj$mle.calc) {
    # if nondiag W is supplied then use it.
    # otherwise assume a diagonal set of weights.
    #
    # NOTE: calculation of  shat involves full set of obs
    # not those colllapsed to the mean.
    if (obj$nondiag.W) {
        shat.GCV <- sqrt(sum((obj$W2 %d*% obj$residuals)^2)/(length(obj$y) - 
            obj$eff.df))
    }
    else {
        shat.GCV <- sqrt(sum((obj$weights * obj$residuals^2)/(length(obj$y) - 
            obj$eff.df)))
    }
    if (mle.calc) {
        rho.MLE <- sum(c(obj$c) * c(obj$yM))/obj$N
        # set rho estimate to zero if negtive. Typically this
        # is an issue of machine precision and very small negative value.
        rho.MLE <- ifelse(rho.MLE < 0, 0, rho.MLE)
        
        #    commented out code for debugging ...
        #      if( rho.MLE< 0) {
        #        stop('problems computing rho.MLE')}
        # commented out is the REML estimate -- lose null space df because of
        # the restiction to orthogonal subspace of T.
        # rhohat<- rho.MLE <- sum(obj$c * obj$yM)/(obj$N - obj$nt)
        # .
        rhohat <- rho.MLE
        shat.MLE <- sqrt(rho.MLE * obj$lambda)
    }
    else {
        rhohat <- rho.MLE <- shat.MLE <- NA
    }
    list(shat.GCV = shat.GCV, rho.MLE = rho.MLE, shat.MLE = shat.MLE, 
        rhohat = rhohat)
}

"Krig.replicates" <- function(out=NULL, x,y, Z=NULL, weights=rep( 1, length(y)),
                               verbose = FALSE) {
    if( is.null(out)){
      out<- list( x=x, y=y, N= length(y), Z=Z, weights=weights)
    }
    rep.info <- cat.matrix(out$x)
    if (verbose) {
        cat("replication info", fill = TRUE)
        print(rep.info)
    }
    # If no replicates are found then reset output list to reflect this condition
    uniquerows <- !duplicated(rep.info)
    if (sum(uniquerows) == out$N) {
        shat.rep <- NA
        shat.pure.error <- NA
        pure.ss <- 0
        # coerce 'y' data vector as a single column matrix
        yM <- as.matrix(out$y)
        weightsM <- out$weights
        xM <- as.matrix(out$x[uniquerows, ])
        # coerce ZM to matrix
        if (!is.null(out$Z)) {
            ZM <- as.matrix(out$Z)
        }
        else {
            ZM <- NULL
        }
    }
    # collapse over spatial replicates
    else {
        rep.info.aov <- fast.1way(rep.info, out$y, out$weights)
        shat.pure.error <- sqrt(rep.info.aov$MSE)
        shat.rep <- shat.pure.error
        # copy  replicate means as a single column matrix
        yM <- as.matrix(rep.info.aov$means)
        weightsM <- rep.info.aov$w.means
        xM <- as.matrix(out$x[uniquerows, ])
        # choose some Z's for replicate group means
        if (!is.null(out$Z)) {
            ZM <- as.matrix(out$Z[uniquerows, ])
        }
        else {
            ZM <- NULL
        }
        pure.ss <- rep.info.aov$SSE
        if (verbose) 
            print(rep.info.aov)
    }
    return(list(yM = yM, xM = xM, ZM = ZM, weightsM = weightsM, 
        uniquerows = uniquerows, shat.rep = shat.rep, shat.pure.error = shat.pure.error, 
        pure.ss = pure.ss, rep.info = rep.info))
}

Krig.transform.xY <- function(obj, knots, verbose = FALSE) {
    # find all replcates and  collapse to unique locations and mean response
    # and pooled variances and weights.
    out <- Krig.replicates(obj, verbose = verbose)
    if (verbose) {
        cat("yM from Krig.transform.xY", fill = TRUE)
        print(out$yM)
    }
    #
    # save information about knots.
    if (is.na(knots[1])) {
        out$knots <- out$xM
        out$mle.calc <- TRUE
        out$knot.model <- FALSE
    }
    else {
        out$mle.calc <- FALSE
        out$knot.model <- TRUE
        out$knots <- knots
    }
    #
    # scale x, knot locations and  save transformation info
    #
    out$xM <- transformx(out$xM, obj$scale.type, obj$x.center, 
        obj$x.scale)
    out$transform <- attributes(out$xM)
    out$knots <- scale(out$knots, center = out$transform$x.center, 
        scale = out$transform$x.scale)
    #
    #
    #verbose block
    #
    if (verbose) {
        cat("transform", fill = TRUE)
        print(out$transform)
    }
    if (verbose) {
        cat("knots in transformed scale", fill = TRUE)
        print(knots)
    }
    return(out)
}

"Krig.updateY" <- function(out, Y, verbose = FALSE, 
    yM = NA) {
    #given new Y values but keeping everything else the same finds the
    #new u vector and pure error SS associated with the Kriging estimate
    # the steps are
    # 1) standardize if neccesary
    # 2) find means, in the case of replicates
    # 3) based on the decomposition, multiply a weighted version of yM
    #    with a large matrix extracted from teh Krig object out.
    #
    # The out object may be large. This function is written so that out is # #not changed with the hope that it is not copied locally  in this  #function .
    # All of the output is accumulated in the list out2
    #STEP 1
    #
    # transform Y by mean and sd if needed
    #
    if (out$correlation.model) {
        Y <- (Y - predict(out$mean.obj, out$x))/predict(out$sd.obj, 
            out$x)
        if (verbose) 
            print(Y)
    }
    #
    #STEP 2
    if (is.na(yM[1])) {
        out2 <- Krig.ynew(out, Y)
    }
    else {
        out2 <- list(yM = yM, shat.rep = NA, shat.pure.error = NA, 
            pure.ss = NA)
    }
    if (verbose) {
        print(out2)
    }
    #
    #STEP3
    #
    # Note how matrices are grabbed from the Krig object
    #
    if (verbose) 
        cat("Type of decomposition", out$decomp, fill = TRUE)
    if (out$decomp == "DR") {
        #
        #
        u <- t(out$matrices$G) %*% t(out$matrices$X) %*% (out$weightsM * 
            out2$yM)
        #
        # find the pure error sums of sqaures.
        #
        temp <- out$matrices$X %*% out$matrices$G %*% u
        temp <- sum((out$W2 %d*% (out2$yM - temp))^2)
        out2$pure.ss <- temp + out2$pure.ss
        if (verbose) {
            cat("pure.ss", fill = TRUE)
            print(temp)
            print(out2$pure.ss)
        }
    }
    #####
    ##### end DR decomposition block
    #####
    ####
    #### begin WBW decomposition block
    ####
    if (out$decomp == "WBW") {
        #### decomposition of Q2TKQ2
        u <- c(rep(0, out$nt), t(out$matrices$V) %*% qr.q2ty(out$matrices$qr.T, 
            out$W2 %d*% out2$yM))
        if (verbose) 
            cat("u", u, fill = TRUE)
        #
        # pure error in this case from 1way ANOVA
        #
        if (verbose) {
            cat("pure.ss", fill = TRUE)
            print(out2$pure.ss)
        }
    }
    #####
    ##### end WBW block
    #####
    out2$u <- u
    out2
}
Krig.which.lambda <- function(out) {
    #
    # determine the method for finding lambda
    #  Note order
    # default is to do 'gcv/REML'
    out2 <- list()
    # copy all all parameters to out2 just to make this
    # easier to read.
    out2$method <- out$method
    out2$lambda.est <- NA
    out2$lambda <- out$lambda
    out2$eff.df <- out$eff.df
    out2$rho <- out$rho
    out2$sigma2 <- out$sigma2
    if (!is.na(out2$lambda) | !is.na(out2$eff.df)) {
        #
        # this indicates lambda has been supplied and leads to
        # the cholesky type computational approaches
        #        -- but only if GCV is FALSE
        #
        out2$method <- "user"
    }
    out2$GCV <- out$GCV
    if (!is.na(out2$eff.df)) {
        #
        # this indicates df has been supplied and needs
        # GCV to be true to compute the lambda
        # that matches the df
        #
        out2$GCV <- TRUE
    }
    if (!is.na(out2$rho) & !is.na(out2$sigma2)) {
        out2$method <- "user"
        out2$lambda <- out2$sigma2/out2$rho
    }
    #
    # NOTE: method='user' means that a value of lambda has been supplied
    #        and so GCV etc to determine lambda is not needed.
    #  gcv TRUE means that the decompositions will be done to
    #    evaluate the estimate at arbitrary lambda (and also be
    #    able to compute the effective degrees of freedom).
    #
    #    The fixed lambda calculations are very efficient but
    #    do not make it feasible for GCV/REML  or effective degrees of
    #    freedom calculations.
    #
    out2$fixed.model <- (out2$method == "user") & (!out2$GCV)
    #
    return(out2)
}

"Krig.ynew" <- function(out, y = NULL, yM = NULL) {
    #
    # calculates the collapsed y (weighted) mean vector based on the
    # X matrix and weights from the out object.
    # or just passes through the collapsed mean data if passed.
    #
    #
    # If there are no replicated obs. then return the full vector
    # pure error ss is zero
    #
    shat.rep <- NA
    shat.pure.error <- NA
    pure.ss <- 0
    # if no y's are given then it is assumed that one should use the
    # yM from the original data used to create the Krig object
    if (is.null(yM) & is.null(y)) {
        yM <- out$yM
    }
    #
    # case when yM is passed no calculations are needed
    #
    if (!is.null(yM)) {
        return(list(yM = as.matrix(yM), shat.rep = NA, shat.pure.error = NA, 
            pure.ss = 0))
    }
    #
    # no reps case
    #
    if (length(unique(out$rep.info)) == out$N) {
        return(list(yM = as.matrix(y), shat.rep = NA, shat.pure.error = NA, 
            pure.ss = 0))
    }
    #
    #  check that y is the right length
    #
    if (length(y) != out$N) {
        stop(" the new y vector is the wrong length!")
    }
    #
    # case when full y data is passed and replicate means need to be found
    #
    if (length(unique(out$rep.info)) < out$N) {
        #
        # calculate means by pooling Replicated obseravations but use the
        # the right weighting.
        #
        rep.info.aov <- fast.1way(out$rep.info, y, out$weights)[c("means", 
            "MSE", "SSE")]
        shat.pure.error <- sqrt(rep.info.aov$MSE)
        shat.rep <- shat.pure.error
        return(list(yM = rep.info.aov$means, shat.rep = shat.rep, 
            shat.pure.error = shat.pure.error, pure.ss = rep.info.aov$SSE))
    }
}
