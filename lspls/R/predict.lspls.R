### predict.lspls: predict method
### $Id: predict.lspls.R 38 2009-07-18 12:19:18Z bhm $

## The plan:  Build up a new `predictors' by calculateing new
## scores, and use object$coefficients to get new predictions.

predict.lspls <- function(object, newdata, type = c("response", "scores"),
                          na.action = na.pass,...) {
    ## Check args:
    type <- match.arg(type)
    ## Build new data matrices:
    tt <- delete.response(terms(object))
    if (missing(newdata) || is.null(newdata))
        mf <- model.frame(object)[-1]
    else {
        mf <- model.frame(tt, newdata, na.action = na.action)
        if (!is.null(cl <- attr(tt, "dataClasses")))
            .checkMFClasses(cl, mf)
    }
    matrices <- apply(attr(tt, "factors"), 2, function(x) mf[,which(x > 0)])
    newX <- as.matrix(matrices[[1]])
    newZ <- matrices[-1]

    ## Parametres:
    nObs <- nrow(newX)

    ## Containers:
    V <- matrix(nrow = nObs, ncol = ncol(object$predictors))

    ## Start with X:
    nVar <- ncol(newX)
    V[,1:nVar] <- newX

    ## Walk through the plsr models:
    for (i in seq(along = object$models)) {
        ##cat("i =", i, "\n")
        M <- newZ[[i]]
        if (is.matrix(M)) {             # Single matrix
            Mo <- M  - V[,1:nVar] %*% object$orthCoefs[[i]]  # Orth. M
            V[,nVar + (1:object$ncomp[[i]])] <-
                sweep(Mo, 2, object$models[[i]]$Xmeans) %*%
                    object$models[[i]]$projection
            nVar <- nVar + object$ncomp[[i]]
        } else {                        # Parallell matrices
            ## The variables to be added in the present step:
            Vadd <- matrix(nrow = nObs, ncol = sum(object$ncomp[[i]]))
            added <- 0
            for (j in seq(along = M)) {
                ##cat("j =", j, "\n")
                ## Walk through Z[[i]]
                Mo <- M[[j]]  - V[,1:nVar] %*% object$orthCoefs[[i]][[j]]
                Vadd[,added + (1:object$ncomp[[i]][[j]])] <-
                    sweep(Mo, 2, object$models[[i]][[j]]$Xmeans) %*%
                        object$models[[i]][[j]]$projection
                added <- added + object$ncomp[[i]][[j]]
            }
            V[,nVar + (1:sum(object$ncomp[[i]]))] <- Vadd
            nVar <- nVar + sum(object$ncomp[[i]])
        } # if
    } # for
    ## Now V contains the new values of the prediction variables
    if (type == "scores")
        return(V)
    else
        return(V %*% object$coefficients)
} # function
