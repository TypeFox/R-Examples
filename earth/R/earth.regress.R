# earth.regress.R: earth.regress is used only for testing Regress in earth.c

earth.regress <- function(
    x         = stop("no 'x' argument"), # NAs are not allowed in x or y
    y         = stop("no 'y' argument"),
    weights   = NULL,               # case weights
    used.cols = NULL)
{
    # following copied from header of earth.fit
    # expand factors, convert to double matrix with col names
    env <- parent.frame()
    x <- expand.arg(x, env, trace=0)
    y <- expand.arg(y, env, trace=0, is.y.arg=TRUE, trunc.deparse(substitute(y)))
    if(nrow(x) == 0)
        stop0("no 'x' values")
    if(ncol(x) == 0)    # this happens for example for earth(Volume~Volume,data=trees)
        stop0("no 'x'")
    if(nrow(x) != nrow(y))
        stop0("nrow(x) ", nrow(x), " != nrow(y) ", nrow(y))
    if(!all(is.double(x)))
        stop0("non double entries in 'x' argument")
    if(!all(is.double(y)))
        stop0("non double entries in 'y' argument")
    if(!is.null(weights) && length(weights) != nrow(x))
        stop0("length(weights) ", length(weights), " != nrow(x) ", nrow(y))

    # add intercept to x
    colnames. <- colnames(x)
    x <- cbind(repl(1, nrow(x)), x)
    colnames(x) <- c("(Intercept)", colnames.)

    nresp <- ncol(y)
    ncols <- ncol(x)
    ncases <- nrow(x)

    weights <- get.weights(weights, ncases)
    if(is.null(used.cols)) {
        used.cols <- repl(TRUE, ncols)
        coefficients <- matrix(1.0, nrow=ncol(x), ncol=nresp)
    } else {
        if(!is.logical(used.cols))
            stop0("used.cols is not logical")
        if(length(used.cols) != ncol(x)-1)     # -1 for intercept added above
            stop0("length(used.cols) != ncol(x)")
        check.index(used.cols, "used.cols", x, is.col.index=TRUE)
        used.cols <- c(TRUE, used.cols)         # add intercept
        coefficients <- matrix(1.0, nrow=ncol(x) - sum(!used.cols), ncol=nresp)
    }
    rownames(coefficients) <- colnames(x)[used.cols]
    colnames(coefficients) <- colnames(y)

    rval <- .C("RegressR",
        coefficients = coefficients, # double  Betas[]     out: nUsedCols * nResp
        residuals = matrix(1.0, nrow=ncases, ncol=nresp),
                                # double       Residuals[] out: nCases * nResp
        rss = double(1),        # double*      pRss        out: RSS, summed over all nResp
        diags = double(ncols),  # double       Diags[]     out:
        rank = integer(1),      # int*         pnRank      out: nbr of indep cols in x
        pivots = integer(ncols),# int          iPivots[]   out: nCols
        as.double(x),           # const double x[]         in: nCases x nCols
        as.double(y),           # const double y[]         in: nCases x nResp
        as.integer(ncases),     # const int*   pnCases     in: number of rows in x and in y
        as.integer(nresp),      # const int*   pnResp      in: number of cols in y
        as.integer(ncols),      # int*         pnCols      in: number of columns in x
        as.integer(used.cols),  # const bool   UsedCols[]) in: specifies used columns in x
        PACKAGE="earth")

    rval$fitted.values <- y - rval$residuals
    rval$call <- match.call()
    rval
}
