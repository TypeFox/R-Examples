"grouped" <-
function(formula, link = c("identity", "log", "logit"), distribution = c("normal", "t", "logistic"), 
                        data, subset, na.action, str.values, df = NULL, iter = 3, ...){
    cl <- match.call()
    mf <- match.call(expand.dots=FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    if(is.empty.model(mt)){
        fit <- list(coefficients = NULL, logLik = NULL, hessian = NULL, k = NULL, n = n)
    } else{
        distr <- match.arg(distribution)
        if(is.null(df) && distr == "t")
            stop("You must specify the degrees of freedom for the Student's-t distribution.\n")
        if(!is.null(df) && distr != "t")
            warning("you specified the `df' argument and you don't use the Student's-t as reference distribution.\n")
        link <- match.arg(link)
        X <- model.matrix(mt, mf)
        fit <- grouped.fit(y, X, link, distr, df, starts = str.values, iter)
    }
    fit$call <- cl
    class(fit) <- "grouped"
    fit
}

