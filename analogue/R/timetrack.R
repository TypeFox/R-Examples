## produce an object that contains an ordination
## and predict new locations for core samples
`timetrack` <- function(X, passive, env,
                        method = c("cca", "rda"),
                        transform = "none",
                        formula, ##type = c("wa","lc"),
                        scaling = 3, rank = "full",
                        ##model = c("CCA", "CA"),
                        join = "left",
                        ...) {
    origX <- X ## store for later
    namX <- deparse(substitute(X))
    namP <- deparse(substitute(passive))
    ## Apply a transformation - let tran deal with arg matching
    if(!isTRUE(all.equal(transform, "none"))) {
        X <- tran(X, method = transform, ...)
        passive <- tran(passive, method = transform, ...)
    }
    ## merge X and passive
    dat <- join(X, passive, type = join)
    X <- dat[[1]]
    passive <- dat[[2]]
    ## common set of species
    tmp <- colSums(X > 0) > 0 ##& colSums(passive > 0) > 0
    X <- X[, tmp]
    passive <- passive[, tmp]
    ## check what type of ordination is required
    ##if(isTRUE(missing(method)))
    ##    method <- "cca"
    method <- match.arg(method)
    FUN <- match.fun(method)
    ## if no env do unconstrained
    if(isTRUE(missing(env))) {
        namE <- NA
        formula <- FALSE
        ord <- FUN(X = X, ...)
    } else {
        namE <- deparse(substitute(env))
        ## check env is same length as nrow(X)
        if(!isTRUE(all.equal(NROW(env), nrow(X))))
            stop("'X' and 'env' imply different numbers of observations")
        ## check if a formula is present
        if(isTRUE(missing(formula))) {
            formula <- FALSE
            ord <- FUN(X = X, Y = env, ...)
        } else {
          mf <- match.call(expand.dots = FALSE)
          m <- match(c("X","passive","env","transform","scaling","rank",
                       "model","method","condition"), names(mf), 0L)
          mf <- mf[-m]
          mf[[1]] <- as.name("model.frame")
          mt <- terms(formula, data = env, simplify = TRUE)
          mf$formula <- formula(mt, data = env)
          mf$data <- env
          dots <- list(...)
          mf[names(dots)] <- NULL
          mf <- eval(mf, parent.frame())
          ord <- FUN(X = X, Y = mf, ...)
        }
    }
    ## process predict args
    ##if(isTRUE(missing(model)))
    ##    model <- "CCA"
    ##model <- match.arg(model)
    ## fitted values for passive
    pred <- predict(ord, newdata = passive, type = "wa",
                    scaling = scaling, model = "CCA", rank = rank)
    pred2 <- predict(ord, newdata = passive, type = "wa",
                     scaling = scaling, model = "CA", rank = rank)
    pred <- cbind(pred, pred2)
    nams <- list(X = namX, passive = namP, env = namE)
    ## return object
    res <- list(ordination = ord, fitted.values = pred,
                method = method, formula = formula, #type = type,
                scaling = scaling, rank = rank, ##model = model,
                labels = nams, call = match.call(),
                X = origX, transform = transform)
    class(res) <- "timetrack"
    res
}

`print.timetrack` <- function(x, ...) {
    cat("\n")
    writeLines(strwrap("Timetrack Ordination", prefix = "\t"))
    cat("\n")
    writeLines(strwrap(pasteCall(x$call)))
    cat("\n")
    writeLines(strwrap("Ordination Output:"))
    ##cat("\n")
    print(x$ordination, ...)
    invisible(x)
}

## require(analogue)
## data(rlgh, swapdiat)
## mod <- timetrack(swapdiat, rlgh, transform = "hellinger",
##                  method = "rda")
## mod
## plot(mod)

