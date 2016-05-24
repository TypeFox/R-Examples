### lspls.R: The user interface function for fitting models
### $Id: lspls.R 39 2009-07-18 12:22:53Z bhm $

lspls <- function(formula, ncomp, data, subset, na.action,
                  model = TRUE, ...)
{
    ## Get the terms
    mt <- terms(formula, keep.order = TRUE)
    ## Get the model frame
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]                   # Retain only the named arguments
    mf[[2]] <- mt                       # Use the terms instead of the
                                        # formula, to keep the ordering.
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    ## Get the data matrices
    Y <- model.response(mf, "numeric")
    if (is.matrix(Y)) {
        if (is.null(colnames(Y)))
            colnames(Y) <- paste("Y", 1:dim(Y)[2], sep = "")
    } else {
        Y <- as.matrix(Y)
        colnames(Y) <- deparse(formula[[2]])
    }
    ## All the predictor matrices, in correct order:
    matrices <- apply(attr(mt, "factors"), 2, function(x) mf[,which(x > 0)])
    X <- as.matrix(matrices[[1]])
    Z <- matrices[-1]
    ## Make sure ncomp is a list, and repeat it as needed:
    ncomp <- rep(as.list(ncomp), length = length(Z))

    ## Fit the model:
    z <- orthlspls.fit(Y, X, Z, ncomp, ...)
    ## Build and return the object:
    class(z) <- "lspls"
    z$fitted.values <- Y - z$residuals
    z$na.action <- attr(mf, "na.action")
    z$ncomp <- ncomp
    z$call <- match.call()
    z$terms <- mt
    if (model) z$model <- mf
    z
}
