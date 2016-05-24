## CA MARCHE !!!
## m111dii <- mhurdle(comics ~ gender + educ  |  log(incum) +
##                    I(log(incum)^2) + I(log(incum)^3) + size | age|log(incum),
##                    data = Comics, corr = "h1", dist = "n", method = 'bfgs', print.level=3)


 ## m111dii <- mhurdle(comics ~ gender + educ  |  log(incum) +
 ##                   I(log(incum)^2) + I(log(incum)^3) + size | age|log(incum),
 ##                   data = Comics, corr = "h1", dist = "b", method = 'bfgs', print.level=3)

mhurdle <- function(formula, data, subset, weights, na.action,
                    start = NULL, dist = c("ln", "tn", "n", "bc", "ihs"), corr = NULL, ...){
    dots <- list(...)
    oldoptions <- options(warn = -1)
    on.exit(options(oldoptions))
    cl <- match.call()
    posT <- as.list(cl) == "T"
    posF <- as.list(cl) == "F"
    cl[posT] <- TRUE
    cl[posF] <- FALSE
    cl.save <- cl
    dist <- match.arg(dist)
    
    ##
    # compute the model.frame and the model.matrix
    ##

    if (!inherits(formula, "Formula")) formula <- Formula(formula)
    if (length(formula)[2] > 4) stop("at most 4 rhs should be provided in the formula")
    mf <- match.call(expand.dots = FALSE) 
    m <- match(c("formula", "data", "subset", "na.action", "weights"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf$formula <- formula
    mf <- eval(mf, parent.frame())
    X1 <- model.matrix(formula, data = mf , rhs = 1)
    X2 <- model.matrix(formula, data = mf , rhs = 2)
    X3 <- model.matrix(formula, data = mf , rhs = 3)
    
    y <- model.response(mf)
    n <- length(y)
    if (length(X1) == 0) X1 <- NULL
    if (length(X3) == 0) X3 <- NULL
    if (length(X2) == 0) stop("the second hurdle (consumption equation) is mandatory")
    h1 <- !is.null(X1)
    h3 <- !is.null(X3)
    
    if (length(formula)[2] == 4){
        X4 <- model.matrix(formula, data = mf, rhs = 4)
        if (length(X4) == 0) X4 <- NULL
    }
    else X4 <- NULL
    
    if (!is.null(corr)){
        if ((h1 + h3) == 2){
            if (nchar(corr) != 3) stop("corr should contain three characters")
            if (!(corr %in% c("dii", "iid"))) stop("corr should currently be either dii or iid")
            if (corr == "dii") corr <- "h1"
            if (corr == "iid") corr <- "h3"
        }
        else{
            if (nchar(corr) != 1) stop("corr should contain just one character")
            if (!(corr %in% c("i", "d"))) stop("corr should be one of d or i")
            if (corr == "i") corr <- NULL
            else{
                if (h1) corr <- "h1"
                if (h3) corr <- "h3"
            }
        }
    }
    # compute the "naive" model
    Pnull <- mean(y == 0)
    if (dist != "ln"){
        Ec <- mean(y[y > 0])
        Vc <- var(y[y > 0])}
    else{
        Ec <- mean(log(y[y > 0]))
        Vc <- var(log(y[y > 0]))
    }
    start.naive <- c(rep(0.1, 1 + h1 + h3), 1)
    moments <- c(Pnull, Ec, Vc)
    dist.naive <- dist
    if (dist %in% c("bc", "ihs")) dist.naive <- "n"
    naive <- maxLik(lnl.naive, start = start.naive,
                    dist = dist.naive, moments = moments,
                    h1 = h1, h3 = h3);
    coef.naive <- naive$est
    logLik.naive <- structure(naive$max * n, nobs = length(y), df = length(coef.naive), class = "logLik")
    naive <- list(coefficients = coef.naive, logLik = logLik.naive, code = naive$code)
    
    #  Special cases where the models can be estimated directly, without
    #  relevant starting values
  
    # model (1l 1t), no censure, estimate the log-linear or truncated
    # model
    if (!h1 && !h3 && dist != "n"){
        if (dist == "ln") result <- lm( log(y) ~ X2 - 1)
        if (dist == "tn") result <- truncreg(y ~ X2 - 1)
        return(result)
    }
    # models (2li, 2ti) h1 without corr, can be estimated simply in two
    # parts using fit.simple.mhurdle() used as starting values for (3d,
    # 4d)
    if (h1 && !h3 &&  ! (dist %in% c("n", "bc", "ihs")) && is.null(corr)){
        result <- fit.simple.mhurdle(X1, X2, y, dist = dist)
        result$naive <- naive
        result$call <- cl.save
        result$model <- mf
        result$formula <- formula
        return(result)
    }
    # Compute the starting values if not provided
    dist.start <- dist
    if (dist %in% c("bc", "ish")) dist.start <- "n"
    if (is.null(start)) start <- start.mhurdle(X1, X2, X3, y, dist.start)
    # in case of heteroscedasctic model, add K4 zeros to the start
    # vector and the intercept should be ln(sigma_o) (not sigma_o)
    # beacause of the exp form
    if (!is.null(X4)){
        sd.int.pos <- ifelse(h1, ncol(X1), 0) + ncol(X2) + ifelse(h3, ncol(X3), 0) + 1
        sd.last.pos <- sd.int.pos - 1 + ncol(X4)
        start[sd.int.pos] <- log(start[sd.int.pos])
        start <- c(start[1:sd.int.pos], rep(0, ncol(X4) - 1))
    }
    if (!is.null(corr)) start <- c(start, 0.1)
    if (dist == "bc") start <- c(start, 1)
    if (dist == "ihs") start <- c(start, .2)
    # Fit the model
    result <- mhurdle.fit(start, X1, X2, X3, X4, y,
                          gradient = TRUE, fit = FALSE,
                          dist = dist, corr = corr, ...)
    result$naive <- naive
    result$call <- cl.save
    result$formula <- formula
    names(result$coefficients) <- colnames(result$vcov) <-
        rownames(result$vcov) <- nm.mhurdle(result)
    result$model <- mf
    result
}

mhurdle.fit <- function(start, X1, X2, X3, X4, y, gradient = FALSE, fit = FALSE,
                        dist = c("ln", "n", "tn", "bc", "ihs"), corr = NULL, ...){
    start.time <- proc.time()
    f <- function(param) mhurdle.lnl(param, X1 = X1, X2 = X2, X3 = X3, X4 = X4, y = y,
                                     gradient = TRUE, fitted = FALSE,
                                     dist = dist, corr = corr)
    check.gradient <- FALSE
    if (check.gradient){
        ngrad <- c()
        oparam <- start
        fo <- f(start)
        agrad <- apply(attr(fo, "gradient"), 2, sum)
        eps <- 1E-05
        for (i in 1:length(start)){
            oparam[i] <- oparam[i] + eps
            ngrad <- c(ngrad, sum((as.numeric(f(oparam)) - fo) / eps))
            oparam <- start
        }
        print(cbind(start, agrad, ngrad))
    }
    maxl <- maxLik(f, start = start,...)
    coefficients <- maxl$estimate
    fitted <- attr(mhurdle.lnl(coefficients, X1 = X1, X2 = X2, X3 = X3, X4 = X4, y = y,
                               gradient = FALSE, fitted = TRUE,
                               dist = dist, corr = corr), "fitted")
    # La ligne ci-dessous renvoie la contribution de chaque obs a la
    # vraisemblance et au gradient (en attribut)
    logLik <- f(coefficients)
    gradi <- attr(logLik, "gradi")
    logLik <- structure(as.numeric(logLik), df = length(coefficients),
                        nobs = length(y), class = "logLik")
    hessian <- maxl$hessian
    # A calculer, les fitted values pour Box-Cox
    convergence.OK <- maxl$code <= 2
    elaps.time <- proc.time() - start.time
    nb.iter <- maxl$iterations
    eps <- with(maxl, gradient %*% solve(- hessian) %*% gradient)
    est.stat <- list(elaps.time = elaps.time,
                     nb.iter = nb.iter,
                     eps = eps,
                     method = maxl$type,
                     message = maxl$message
                     )
    class(est.stat) <- "est.stat"
    if (! is.null(X4)) sd.names <- colnames(X4) else sd.names <- "sd"
    if (!is.null(corr)) rho.names <- ifelse(corr == "h1", "corr12", "corr13") else rho.names <- NULL
    if (dist %in% c("bc", "ihs")) tr.names <- "tr" else tr.names <- NULL
    coef.names <- list(h1   = colnames(X1),
                       h2   = colnames(X2),
                       h3   = colnames(X3),
                       sd   = sd.names,
                       corr = rho.names,
                       tr   = tr.names)
    result <- list(coefficients  = coefficients,
                   vcov          = - solve(maxl$hessian),
                   fitted.values = fitted,
                   logLik        = logLik,
                   gradient      = gradi,
                   formula       = NULL,
                   model         = NULL,
                   coef.names    = coef.names,
                   call          = NULL,
                   est.stat      = est.stat,
                   naive         = NULL
                   )
    if (ncol(X2) > 1) class(result) <- c("mhurdle")
    result
}
 
# Compute the starting values of the model
start.mhurdle <- function(X1, X2, X3, y, dist){
    h1 <- !is.null(X1)
    h3 <- !is.null(X3)
    # for models (2ld, 2td), estimates of models (2li, 2ti) which can be
    # estimated simply in two parts using fit.simple.mhurdle() are used
    # as starting values
    if (h1 && !h3 &&  dist != "n"){
        result <- fit.simple.mhurdle(X1, X2, y, dist = dist)
        start <- result$coefficients
    }
    # model (3) tobit : use linear model as starting values
    if (!h1 && !h3 &&  dist == "n"){
        lin <- lm(y ~ X2 - 1)
        start <- c(coef(lin), summary(lin)$sigma)
    }
    # model (4, 7) h3 whithout h1
    if (!h1 && h3){
        probit <- glm( (y != 0) ~ X3 - 1, family = binomial(link = 'probit'))
        bX3 <- as.numeric(crossprod(coef(probit), t(X3)))
        Phi3 <- pnorm(bX3)
        yPP <- y * Phi3
        lin <- switch(dist,
                      "ln" = lm(log(yPP)  ~ X2 - 1, subset = y != 0),
                      "n" = lm(yPP       ~ X2 - 1, subset = y != 0),
                      "tn" = truncreg(yPP ~ X2 - 1, subset = y != 0)
                      )
        if (dist %in% c("n", "ln")) start <- c(coef(lin), coef(probit), summary(lin)$sigma)
        if (dist == "tn") start <- c(coef(lin)[- length(coef(lin))], coef(probit), coef(lin)[length(coef(lin))])
    }
    # model (5), double hurdle use model (3i) as starting values
    if (h1 && !h3 && dist == "n"){
        result <- fit.simple.mhurdle(X1, X2, y, dist = dist)
        start <- result$coefficients
    }
    # model (6 and 8)
    if (h1 && h3){
        probit.h3 <- glm( (y != 0) ~ X3 - 1, family = binomial(link = 'probit'))
        probit.h1 <- glm( (y != 0) ~ X1 - 1, family = binomial(link = 'probit'))
        beta3 <- coef(probit.h3)
        beta1 <- coef(probit.h1)
        bX3 <- as.numeric(crossprod(beta3, t(X3)))
        Phi3 <- pnorm(bX3)
        bX1 <- as.numeric(crossprod(beta1, t(X1)))
        P0 <- mean(y > 0)
        yPP <- y * Phi3
        lin <- switch(dist,
                      "ln" = lm(log(yPP)  ~ X2 - 1, subset = y!= 0),
                      "n" = lm(yPP       ~ X2 - 1, subset = y!= 0),
                      "tn" = truncreg(yPP ~ X2 - 1, subset = y!= 0)
                      )
        if(dist != "tn")
            start <- c(beta1, coef(lin), beta3, summary(lin)$sigma)
        else
            start <- c(beta1, coef(lin)[- length(coef(lin))],
                       beta3, coef(lin)[length(coef(lin))])
    }
    start
}

