hdlm <-
function(formula, data, subset, bootstrap = 10, siglevel = 0.05,
    alpha = 0.5, M = NULL, N = NULL, model = TRUE, x = FALSE,
    y = FALSE, scale=TRUE, pval.method=c('median', 'fdr', 'holm', 'QA'),
    ..., FUNCVFIT = NULL, FUNLM = NULL, bayes=FALSE, bayesIters=NULL,
    bayesTune=c(1,1), refit=FALSE) 
{

    # Check inputs for correct form
    if( !is.numeric(bootstrap) | !is.numeric(siglevel) | !is.numeric(alpha)) {
        stop('variables bootstrap, siglevel, and alpha must all be of type numeric')
    }
    if( round(bootstrap) != bootstrap | bootstrap <= 0 | length(bootstrap) != 1 ) {
        stop('bootstrap must be an integer greater than or equal to 1') 
    }
    if( siglevel <= 0 | siglevel >= 1 | length(siglevel) != 1 ) {
        stop('siglevel must be a number strictly between 0 and 1') 
    }
    if( alpha <= 0 | alpha > 1 | length(alpha) != 1 ) {
        stop('alpha must be a number greater than 0 and no more than 1') 
    }
    if( !is.logical(model) | !is.logical(x) | !is.logical(y) | !is.logical(scale) | !is.logical(bayes)) {
        stop('model, x, y, scale, and bayes must all be of type logical')
    }
    if( !is.logical(refit)) {
        if(!is.numeric(refit) | length(refit) != 1) stop('refit must be a length one vector of class logical or numeric')
        if( refit <= 0 | refit > 1 ) stop('if refit is a numeric, it must be positive and no more than one')
    }
    if(!is.null(M) ) {
        if(!is.numeric(M) | M != round(M) | M < 0 ) {
            warning('M must be non-negative integer or NULL; setting to NULL')
            M <- NULL
        }
    }
    if(!is.null(N) ) {
        if(!is.numeric(N) | N != round(N) | N <= 0 ) {
            warning('N must be positive integer or NULL; setting to NULL')
            N <- NULL
        }
    }
    if( !is.function(FUNLM) & !is.null(FUNLM) ) {
        warning('FUNLM must be function or NULL; setting to NULL')
        FUNLM <- NULL
    }
    if( !is.function(FUNCVFIT) & !is.null(FUNCVFIT) ) {
        warning('FUNCVFIT must be function or NULL; setting to NULL')
        FUNCVFIT <- NULL
    }
    pval.method <- pval.method[[1]]
    if(!(pval.method %in% c('median', 'fdr', 'holm', 'QA'))) {
        warning('Unrecognized pval.method; setting to "mean"')
        pval.method <- 'mean'
    }
    if(pval.method == 'QA' & bootstrap < 20) {
        warning('QA not supported for fewer than 20 bootstrap runs; setting to "mean"')
        pval.method <- 'mean'
    }

    # Parse formula and options into model matrix. Similar to lm() code
    ret.x <- x
    ret.y <- y
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    w <- as.vector(model.weights(mf))
    x <- sparse.model.matrix(mt, mf, contrasts)
    if(!is.matrix(x)) x <- as.matrix(x)

    # Generally, the intercept should not be penalized; therefore the intercept 
    # is taken out of the data matrix 'x' if included; most selectors (including
    # the default glmnet) include the offset as a natural addition. Although, should
    # be left in when using Bayesian method:
    intercept <- FALSE
    if(colnames(x)[[1L]] == "(Intercept)" & !bayes) {
      x <- x[,-1] 
      intercept <- TRUE
    }

    # If M = 0, don't do second stage; just return one model selection output
    p <- ncol(x)
    n <- nrow(x)
    if(is.null(N)) N <- floor(n/2)
    if(is.null(M)) M <- floor((n - N) * 0.9)

    if(M == 0) {
        z <- hdlm.fit(x, y, bootstrap=1, siglevel, intercept, alpha = 0.5, M = 0, N = nrow(x),
                   sd.off = 0,scale=scale, pval.method, ..., FUNCVFIT = FUNCVFIT, FUNLM = FUNLM)
        if(M == 0) return(z)
    }

    # If not bayes, pass data to hdlm.fit; otherwise pass data to bayes.hdlm.fit
    if(!bayes) {
        z <- hdlm.fit(x, y, bootstrap, siglevel, intercept, alpha = 0.5, M = M, N = N,
                   sd.off = 0,scale=scale, pval.method, ..., FUNCVFIT = FUNCVFIT, FUNLM = FUNLM)
    } else {
        if(!is.numeric(bayesTune) | length(bayesTune) != 2 | min(bayesTune) < 0 | max(bayesTune) > 1) {
          warning('Invalid bayesTune parameter; see help pages for details')
        }
        z <- bayes.hdlm.fit(x, y, siglevel=siglevel, bayesIters=bayesIters, bayesTune=bayesTune)
        pval.method <- 'mean'
    }

    # Refit or not:
    if( refit != FALSE) {
        refit <- as.numeric(refit)
        model <- which(z$p.value < refit)
        x <- x[,model]
        if(is.null(FUNLM)) FUNLM <- lm
        if(intercept) {
            return(FUNLM(y ~ x))
        } else {
            return(FUNLM(y ~ x - 1))
        }
    } else {
        # Essentially the same code as lm()
        class(z) <- c("hdlm")
        z$pval.method <- pval.method
        if(intercept & !bayes) {
            names(z$coefficients) <- c('(Intercept)', colnames(x))
        } else {
        names(z$coefficients) <- colnames(x)
        }
        z$xlevels <- .getXlevels(mt, mf)
        z$call <- cl
        z$terms <- mt
        if (model) 
            z$model <- mf
        if (ret.x) 
            z$x <- x
        if (ret.y) 
            z$y <- y
        z$level <- 1
        z$siglevel <- siglevel
        return(z)
    }
}

