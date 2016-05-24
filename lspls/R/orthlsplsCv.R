### orthlsplsCv.R: Cross-validation, orthogonalizing version
###
### $Id: orthlsplsCv.R 35 2009-07-18 11:29:48Z bhm $

## The algorithm is based on recursion, after X has been handled.

orthlsplsCv <- function(Y, X, Z, ncomp, segments, trace = FALSE, ...) {

    ## The recursive function:
    ## It uses the following variables from orthlsplsCv
    ## - Z: list of spectral matrices
    ## - ncomp: list of #comps to use in the CV
    ## - segment: indices of the segment to be predicted
    ## - cvPreds: array of predictions; dim: c(nObs, nResp, unlist(ncomp))
    ## - pls.fit: the pls fit function
    cvPredRest <- function(indices, prevCalib, prevPred, prevComps, prevRes) {
        ## indices is the indices of the remaining matrices
        ## prevCalib is a matrix with the X vars and scores used in the
        ##   previous calibrations
        ## prevPred is a matrix with the X vars and scores used in the
        ##   previous predictions
        ## prevComps is the numbers of components used in the previous
        ##   calibrations
        ## prevRes is the residuals from the previous calibrations

        ## The general idea is to handle the first matrix or list of
        ## matrices (Z[[indices[1]]]), and then recall itself on the rest of
        ## the matrices.

        ## The matrix/matrices to handle in this call:
        ind <- indices[1]
        M <- Z[[ind]]
        if (is.matrix(M)) {             # A single matrix
            ## Orthogonalise the calibration spectra wrt. prevVars
            Mcal <- M[-segment,]
            Mo <- orth(Mcal, prevCalib)

            ## Orthogonalise the prediction spectra
            Mpred <- M[segment,]
            Mpo <- Mpred - prevPred %*% Corth(Mcal, prevCalib)
            ## mal: Zorig[i,] - Xorig[i,] %*% Co(Xorig[-i,]) %*% Zorig[-i,]

            ## Estimate a model prevRes ~ orth. spectra + res
            plsM <- pls.fit(Mo, prevRes, ncomp[[ind]])
            ## Save scores:
            calScores <- plsM$scores

            ## Predict new scores and response values
            predScores <- sweep(Mpo, 2, plsM$Xmeans) %*% plsM$projection
            ## FIXME: Only for orth.scores alg:
            predVals <- array(dim = c(nrow(predScores), dim(plsM$Yloadings)))
            for (a in 1:ncomp[[ind]])
                predVals[,,a] <-
                    sweep(predScores[,1:a] %*% t(plsM$Yloadings[,1:a, drop=FALSE]),
                          2, plsM$Ymeans, "+")

            ## Add the predictions to the outer cvPreds variable
            ## Alt. 1:  Calculate the 1-index indices manuall, and use
            ## single indexing (probably quickest, but requires a loop).
            ## Alt. 2:  Use matrix indexing with an expanded grid.
            ##eg <- expand.grid(segment, 1:nResp, ncomp[[ind]])
            ##indMat <- do.call("cbind", c(eg[1:2], as.list(prevComps), eg[3]))
            ##cvPreds[indMat] <- cvPreds[indMat] + predVals
            ## Alt. 3: Build and eval an expression which does what we want:
            ncomps <- length(prevComps)
            nrest <- length(dim(cvPreds)) - ncomps - 3
            dummy <- Quote(cvPreds[segment,])
            dummy[4 + seq(along = prevComps)] <- prevComps + 1
            dummy[5 + ncomps] <- -1
            if (nrest > 0) dummy[5 + ncomps + 1:nrest] <- dummy[rep(4, nrest)]
            eval(substitute(dummy <<- dummy + c(predVals), list(dummy = dummy)))

            ## Return if this is the last matrix/set of matrices
            if (length(indices) == 1) return()

            ## Calculate new residuals
            newResid <- - plsM$fitted.values + c(prevRes)

            ## To save space: drop the model object(s)
            rm(plsM)

            ## Recursively call ourself for each number of components in the
            ## present model
            for (i in 0:ncomp[[ind]])
                Recall(indices[-1], # Remove the index of the current matrix
                       cbind(prevCalib, calScores[,seq(length = i)]), # Add the scores we've used
                       cbind(prevPred, predScores[,seq(length = i), drop=FALSE]), # Add the scores we've predicted
                       c(prevComps, i), # and the number of comps
                       if (i > 0) newResid[,,i] else prevRes) # update the residual

        } else {                        # List of parallell matrices
            Scal <- list()              # The current calibration scores
            Spred <- list()             # The current prediction scores
            for (j in seq(along = M)) {
                ## Orthogonalise the calibration spectra wrt. prevVars
                Mcal <- M[[j]][-segment,]
                Mo <- orth(Mcal, prevCalib)

                ## Orthogonalise the prediction spectra
                Mpred <- M[[j]][segment,]
                Mpo <- Mpred - prevPred %*% Corth(Mcal, prevCalib)

                ## Estimate a model prevRes ~ orth. spectra + res
                plsM <- pls.fit(Mo, prevRes, ncomp[[ind]][j])
                ## Save scores:
                Scal[[j]] <- plsM$scores

                ## Predict new scores
                Spred[[j]] <- sweep(Mpo, 2, plsM$Xmeans) %*% plsM$projection
            }
            ## To save space: drop the model object
            rm(plsM)

            ## Loop over the different combinations of #comps:
            nComps <- expand.grid(lapply(ncomp[[ind]], seq, from = 0))
            for (cind in 1:nrow(nComps)) {
                newComps <- nComps[cind,]
                comps <- c(prevComps, unlist(newComps))
                ## Predict new response values
                calScores <-
                    do.call("cbind", mapply(function(B, b) B[,seq(length=b), drop=FALSE],
                                            Scal, newComps, SIMPLIFY = FALSE))
                predScores <-
                    do.call("cbind", mapply(function(B, b) B[,seq(length=b), drop=FALSE],
                                            Spred, newComps, SIMPLIFY = FALSE))
                if (all(newComps == 0)) {
                    newResid <- prevRes
                } else {
                    lsS <- lm.fit(calScores, prevRes) # FIXME: How about intercept/numerical accurracy?
                    newResid <- lsS$residuals
                    predVals <- predScores %*% lsS$coefficients
                    rm(lsS)
                    ## Add the predictions to the outer cvPreds variable.  Build
                    ## and eval an expression which does what we want:
                    nc <- length(comps)
                    nrest <- length(dim(cvPreds)) - nc - 2
                    dummy <- Quote(cvPreds[segment,])
                    dummy[4 + seq(along = comps)] <- comps + 1
                    if (nrest > 0) dummy[4 + nc + 1:nrest] <- dummy[rep(4, nrest)]
                    eval(substitute(dummy <<- dummy + c(predVals),
                                    list(dummy = dummy)))
                }

                if (length(indices) > 1) { # There are more matrices to fit
                    ## Recursively call ourself
                    Recall(indices[-1], # Remove the index of the current matrices
                           cbind(prevCalib, calScores), # Add the scores we've used
                           cbind(prevPred, predScores), # Add the scores we've predicted
                           c(comps), # and the number of comps
                           newResid) # use the new residual

                }
            } ## for
        } ## if
    } ## recursive function

    ## Setup:
    nObs <- nrow(X)
    nResp <- ncol(Y)
    ## cvPreds: the cross-validated predictions:
    cvPreds <- array(0, dim = c(nObs, nResp, unlist(ncomp) + 1))
    ## Build an unevaluated expression that will insert the predictions into
    ## cvPreds[segment,,...,] when evaluated:
    ndim <- length(dim(cvPreds))
    ## This creates an expression with the empty index argument repeated as
    ## many times as neccessary:
    dummy <- Quote(cvPreds[segment,])[c(1:3, rep(4, ndim - 1))]
    ## Substitute this in an assignment statement:
    addPredictions <- substitute(dummy <- predVals, list(dummy = dummy))

    ## Choose PLS fit algorithm.  FIXME: Support other algs?
    pls.fit <- oscorespls.fit

    ## The main cross-validation loop
    trace <- isTRUE(trace)
    if (trace) cat("Segment: ")
    for (n.seg in seq_along(segments)) {
        if (trace) cat(n.seg, "")
        segment <- segments[[n.seg]]

        ## Handle X
        lsX <- lm.fit(X[-segment,, drop = FALSE], Y[-segment,, drop = FALSE])
        resid <- lsX$residuals
        predVals <- X[segment,, drop = FALSE] %*% lsX$coefficients

        ## Insert the predictions into the cvPred array:
        eval(addPredictions)

        ## Handle the rest of the matrices:
        cvPredRest(indices = 1:length(ncomp),
                   prevCalib = X[-segment,, drop = FALSE],
                   prevPred = X[segment,, drop = FALSE],
                   prevComps = c(),
                   prevRes = resid)
    }
    if (trace) cat("\n")

    dimnames(cvPreds) <- c(list(1:nObs, colnames(Y)),
                           lapply(unlist(ncomp), function(x) 0:x))
    return(cvPreds)
} ## function
