ltm <-
function (formula, constraint = NULL, IRT.param, start.val = NULL, na.action = NULL, control = list()) {
    cl <- match.call()
    tm <- terms(formula)
    av <- all.vars(formula)
    X <- eval(attr(tm, "variables"), list(z1 = 0, z2 = 0), parent.frame())[[1]]
    oX <- X
    if ((!is.data.frame(X) & !is.matrix(X)) || ncol(X) == 1)
        stop("\nthe left-hand side of 'formula' must be either a numeric matrix or a data.frame, with at least two columns.\n")    
    X <- data.matrix(X)
    if (any(its <- apply(X, 2, function (x) { x <- x[!is.na(x)]; length(unique(x)) } ) > 2))
        stop("'data' contain more that 2 distinct values for item(s): ", paste(which(its), collapse = ", "))
    X <- apply(X, 2, function (x) if (all(unique(x) %in% c(1, 0, NA))) x else x - 1)    
    if (!is.null(na.action))
        X <- na.action(X)
    factors <- sum(av %in% c("z1", "z2"))
    if (factors > 2)
        stop("\nMaximum number of factors to include is 2.")
    if ((factors == 1 & !"z1" %in% av) || (factors == 2 & !c("z1", "z2") %in% av))  
        stop("\nyou have to use 'z1' for the first factor and 'z2' for the second one.")
    tm.lab <- attr(tm, "term.labels")
    ltst <- list(factors = factors, inter = "z1:z2" %in% tm.lab, quad.z1 = "I(z1^2)" %in% tm.lab, 
                 quad.z2 = "I(z2^2)" %in% tm.lab)
    IRT.param <- if (missing(IRT.param)) {
        factors == 1 & !ltst$quad.z1
    } else {
        if (!is.logical(IRT.param))
            stop("'IRT.param' must be logical")
        if (IRT.param && !(factors == 1 & !ltst$quad.z1)) {
            warning("conversion to the IRT parameterization works only for the two-parameter logistic model.\n")
            FALSE
        } else
            IRT.param
    }
    p <- ncol(X)
    q. <- 1 + factors + sum(unlist(ltst[-1]))
    betas <- start.val.ltm(start.val, X, factors, formula)
    if (!is.null(constraint)) {
        if ((!is.numeric(constraint) | !is.matrix(constraint)) || (nrow(constraint) > p * q. - 1 | ncol(constraint) != 3))
            stop("'constraint' should be a 3-column numeric matrix with at most ", p * q. - 1, " rows (read help file).\n")
        if (any(constraint[, 1] < 1 | constraint[, 1] > p))
            stop("the 1st column of 'constraint' denotes the items and it should between 1 and ", p, " (read help file).\n")
        if (any(constraint[, 2] < 1 | constraint[, 2] > q.))
            stop("the 2nd column of 'constraint' denotes either the intercept or the factor loadings and it should between 1 and ", 
                    q., " (read help file).\n")
        constraint <- constraint[order(constraint[, 1]), , drop = FALSE]
        constraint[, 1:2] <- round(constraint[, 1:2])
        betas[constraint[, 1:2]] <- constraint[, 3]
    }
    con <- list(iter.em = 40, iter.qN = 150, GHk = if (factors == 1) 21 else 15, method = "BFGS", parscale = NULL,
                verbose = getOption("verbose"))
    con[names(control)] <- control
    fit <- ltm.fit(X, betas, constraint, formula, con)
    ltst$nams <- colnames(fit$coefficients)
    fit$ltst <- ltst
    fit$X <- oX
    fit$control <- con
    fit$IRT.param <- IRT.param
    fit$constraint <- constraint
    fit$formula <- formula
    fit$call <- cl
    class(fit) <- "ltm"
    fit
}
