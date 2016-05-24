lqa.formula <-
function (formula, data = list (), weights = rep (1, nobs), subset, na.action, start = NULL, etastart, mustart, offset, ...)
{
    call <- match.call ()


### Extract data from formula environment:
### --------------------------------------

    if (missing (data)) 
      data <- environment (formula)

    mf <- match.call (expand.dots = FALSE)
    m <- match (c ("formula", "data", "subset", "weights", "na.action", "etastart", "mustart", "offset"), names (mf), 0)
    mf <- mf[c (1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name ("model.frame")
    mf <- eval (mf, parent.frame())
    mt <- attr (mf, "terms")
    Y <- model.response (mf, "any")

    if (length (dim (Y)) == 1) 
    {
        nm <- rownames (Y)
        dim (Y) <- NULL
        if (!is.null (nm)) 
            names(Y) <- nm
    }
  
    if (length (dim (Y)) > 1)
      stop ("Response y must be univariate")

    X <- if (!is.empty.model (mt))
           model.matrix (mt, mf, contrasts)
         else
           matrix (, nrow (Y), 0)

    intercept <- attr (mt, "intercept") > 0
    nobs <- length (Y)

    if (!is.null (weights) && !is.numeric (weights)) 
        stop("'weights' must be a numeric vector")

    offset <- as.vector (model.offset (mf))

    if (!is.null (weights) && any (weights < 0)) 
        stop("negative weights not allowed")

    if (!is.null(offset)) 
    {
        if (length(offset) == 1) 
           offset <- rep (offset, nrow (Y))
        else 
           if (length (offset) != nrow (Y)) 
            stop (gettextf ("number of offsets is %d should equal %d (number of observations)", 
                length (offset), nrow (Y)), domain = NA)
    }

    mustart <- model.extract (mf, "mustart")
    etastart <- model.extract (mf, "etastart")


### Calling quadpen.fit and preparing the return object:
### ----------------------------------------------------

    fit <- lqa.default (x = X, y = Y, intercept = intercept, weights = weights, start = start, etastart = etastart, mustart = mustart, offset = offset, ...)
    fit$na.action <- attr (mf, "na.action")
    fit <- c (fit, list (call = call, formula = formula, terms = mt, data = data, offset = offset, contrasts = attr (X, "contrasts")))
    class (fit) <- c ("lqa", "glm", "lm")
    fit
}

