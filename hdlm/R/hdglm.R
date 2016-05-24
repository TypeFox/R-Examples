hdglm <-
function(formula, data, subset, family = c("gaussian","binomial","poisson"),
    bootstrap = 10, siglevel = 0.05,
    alpha = 0.5, M = NULL, N = NULL, model = TRUE, x = FALSE,
    y = FALSE, scale=TRUE, pval.method=c('median', 'fdr', 'holm', 'QA'),
    ..., FUNCVFIT = NULL, FUNLM = NULL, bayes=FALSE,
    bayesIters=NULL, bayesTune=NULL, refit=FALSE) 
{

    if(refit != FALSE & bayes == TRUE) {
      warning('hdglm does not support options bayes=TRUE and refit = TRUE; if desired refit manually') 
      refit <- FALSE
    }

    family = family[1]
    familyFUN <- switch(family, gaussian = gaussian,
                            binomial = binomial,
                            poisson = poisson)
    GLMcl <- match.call()

    if(is.null(FUNLM)) {
        FUNLM <- function(formula) return(glm(formula, family=familyFUN))
    }
    if(is.null(FUNCVFIT)) {
        FUNCVFIT <- function(x,y) return(mod.cv.glmnet(x,y,alpha=alpha,standardize=scale,nfolds=3,family=family))
    }

    # If wanting binomial latent Bayes model, run it here without ever passing to hdlm():
    if(bayes & family == 'binomial') {
        # Parse Inputes
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
        x <- model.matrix(mt, mf, contrasts)
        if(is.null(bayesTune)) bayesTune <- 0.1
        if(bayesTune >= 1 | bayesTune <= 0) stop('bayesTune must be strictly between 0 and 1')

        # Make sure Y is 0/1 data
        y <- as.numeric(y)
        if(sum(y == 0 | y == 1) != length(y)) stop('response vector must 0s and 1s for family="binomial"')

        # Call Gibbs sampler:
        z <- bayes.hdglm.fit(x=x,y=y,siglevel=0.05,bayesIters = bayesIters, bayesTune = bayesTune)
    
        # Clean response and return here:
        z$family <- 'binomial'
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
        z$pval.method <- 'bayes' # to print / summarize correctly
        class(z) <- c('hdglm', 'hdlm')
        return(z)
    } else if(bayes & family == 'poisson') {
        # No implemented Bayes solution for Poisson regresssion
        stop('Bayes solution not implemented for poisson regression')
    }
    #  If family == 'Gaussian', and bayes=TRUE, handle via hdlm in normal way

    if(missing(data) & missing(subset)) {
        z <- hdlm(formula=formula,  bootstrap = bootstrap, siglevel = siglevel,
              alpha = alpha, M = M, N = N, model = model, x = x, y = y, scale=scale,
              pval.method=pval.method, ... , FUNCVFIT = FUNCVFIT, FUNLM = FUNLM,
              bayes=bayes, bayesIters=bayesIters, bayesTune=bayesTune, refit=FALSE)
    } else if(missing(subset)) {
        z <- hdlm(formula=formula, data=data, bootstrap = bootstrap, siglevel = siglevel,
              alpha = alpha, M = M, N = N, model = model, x = x, y = y, scale=scale,
              pval.method=pval.method, ... , FUNCVFIT = FUNCVFIT, FUNLM = FUNLM,
              bayes=bayes, bayesIters=bayesIters, bayesTune=bayesTune, refit=FALSE)
    } else {
        z <- hdlm(formula=formula, data=data, subset=subset, bootstrap = bootstrap, siglevel = siglevel,
              alpha = alpha, M = M, N = N, model = model, x = x, y = y, scale=scale,
              pval.method=pval.method, ... , FUNCVFIT = FUNCVFIT, FUNLM = FUNLM,
              bayes=bayes, bayesIters=bayesIters, bayesTune=bayesTune, refit=FALSE)
    }

    if(is.null(M)) M <- 1
    if(M == 0) return(z) 

    # Fix some of the output, as model is not linear if family != 'gaussian':
    z$family <- family
    if(family == 'binomial') {
        vals <- z$fitted - z$resid
        z$sigma.hat <- NULL
        z$fitted <- 1 / (1 + exp(z$fitted))
        z$resid <- z$fitted - vals
    } else if (family == 'poisson') {
        vals <- z$fitted - z$resid
        z$sigma.hat <- NULL
        z$fitted <- exp(z$fitted)
        z$resid <- z$fitted - vals
    }
    z$call <- GLMcl
    if(family != 'gaussian') class(z) <- c('hdglm', class(z))
    return(z)
}

