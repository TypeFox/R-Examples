
# $Id: pqfunctions.R 409 2014-12-01 15:58:09Z thothorn $

pqglht <- function(object) 
{
    betahat <- coef(object)
    covm <- vcov(object)
    m <- coef(object, rhs = TRUE)
    df <- object$df

    ses <- sqrt(diag(covm))
    tstat <- (betahat - m) / ses
    cr <- cov2cor(covm)
    dim <- ncol(cr)

    ### p value function
    pfunction <- function(type = c("univariate", "adjusted", p.adjust.methods), 
                          ...) {

        type <- match.arg(type)

        pfct <- function(q) {
            switch(object$alternative, "two.sided" = {
                      low <- rep(-abs(q), dim)
                      upp <- rep( abs(q), dim)
               }, "less" = {
                      low <- rep(q, dim)
                      upp <- rep(Inf, dim)
               }, "greater" = {
                      low <- rep(-Inf, dim)
                      upp <- rep(q, dim)
               })
               pmvt(lower = low, upper = upp, df = df, corr = cr, ...)
        }

        switch(object$alternative, "two.sided" = {
            if (df > 0) pvals <- 2*(1 - pt(abs(tstat),df))     
            else        pvals <- 2*(1 - pnorm(abs(tstat)))
        }, "less" = {
            if (df > 0) pvals <- pt(tstat,df) 
            else        pvals <- pnorm(tstat)
        }, "greater" = {
            if (df > 0) pvals <- 1 - pt(tstat,df)
            else        pvals <- 1 - pnorm(tstat)
        })

        if (type == "univariate")
            return(pvals)

        if (type == "adjusted") {
            ret <- numeric(length(tstat))
            error <- 0
            for (i in 1:length(tstat)) {
                tmp <- pfct(tstat[i])
                if (attr(tmp, "msg") != "Normal Completion" &&  
                    length(grep("^univariate", attr(tmp, "msg"))) == 0)
                    warning(attr(tmp, "msg"))
                if (error < attr(tmp, "error")) 
                    error <- attr(tmp, "error")
                ret[i] <- tmp
            }
            ret <- 1 - ret
            attr(ret, "error") <- error
            return(ret)
        }

        return(p.adjust(pvals, method = type))
    }

    ### quantile function
    qfunction <- function(conf.level, adjusted = TRUE, ...) {

        tail <- switch(object$alternative, "two.sided" = "both.tails",
                                    "less"      = "lower.tail",
                                    "greater"   = "upper.tail")
        if (adjusted) {
            calpha <- qmvt(conf.level, df = df, corr = cr, tail = tail, 
                           ...)
        } else {
            calpha <- qmvt(conf.level, df = df, corr = matrix(1), tail = tail, ...)
        }
        ret <- calpha$quantile
        attr(ret, "error") <- calpha$estim.prec
        return(ret)
    }
    RET <- list(pfunction = pfunction, qfunction = qfunction,
                coefficients = betahat, sigma = ses, tstat = tstat)
    class(RET) <- "pqglht"
    RET
}

### functions for summary(..., test = ) argument
### univariate p values for each linear hypothesis
univariate <- function() 
{
    function(object) {
        RET <- pqglht(object)
        RET$pvalues <- RET$pfunction("univariate")
        RET$type <- "univariate"
        class(RET) <- "mtest"
        RET
    }
}

### global classical Chisq or F tests
global <- function(type = c("Chisq", "F")) 
{
    type <- match.arg(type)
    
    fct <- function(object) {

        RET <- pqglht(object)
        betahat <- RET$coefficients
        m <- coef(object, rhs = TRUE)
        covm <- vcov(object)

        tmp <- betahat - m
        MP <- MPinv(covm)
        SSH <- t(tmp) %*% MP$MPinv %*% tmp

        q <- MP$rank
        if (type == "F") {
            df <- df.residual(object$model)
            if (is.null(df)) {
                type <- "Chisq"
                warning(sQuote("df.residual"), " is not available for ",
                        sQuote("model"), " a Chisq test is performed ",
                        "instead of the requested F test.")
            }
        }
        if (type == "Chisq") {
            pval <- pchisq(SSH, q, lower.tail = FALSE)
        } else {
            pval <- pf(SSH/q, q, df, lower.tail = FALSE)
        }
        RET$pvalue  <- pval
        RET$type <- type
        RET$SSH <- SSH
        RET$fstat <- SSH/q
        RET$df <- c(q, df)
        class(RET) <- "gtest"
        return(RET)
    }
    return(fct)
}

Ftest <- function() global("F")
Chisqtest <- function() global("Chisq")

### p values adjusted for simultaneous inference
adjusted <- function(type = c("single-step", "Shaffer", "Westfall", "free", p.adjust.methods), 
                     ...) 
{
    type <- match.arg(type)

    ### usual max-type adjustment over all linear hypotheses
    if (type == "single-step") {
        return(function(object) {
            RET <- pqglht(object)
            RET$pvalues <- RET$pfunction("adjusted", ...)
            RET$type <- type
            class(RET) <- "mtest"
            RET
        })
    }

    if (type == "free") {
        return(function(object) {
            K <- object$linfct
            rhs <- object$rhs
            pv <- matrix(0, nrow = nrow(K), ncol = nrow(K))
            colnames(pv) <- rownames(K)
            RET <- pqglht(object)
            for (i in 1:nrow(K)) {
                object$linfct <- K
                object$rhs <- rhs
                tmp <- pqglht(object)
                p <- tmp$pfunction("adjusted", ...)
                pv[i,rownames(K)] <- min(p)
                K <- K[-which.min(p),, drop = FALSE]
                rhs <- rhs[-which.min(p)]
            }
            RET$pvalues <- apply(pv, 2, max)
            RET$type <- type
            class(RET) <- "mtest"
            RET
        })
    }

    ### Westfall (1997, JASA): constraints and correlations 
    ### or
    ### Shaffer (1886, JASA): constraints
    if (type %in% c("Shaffer", "Westfall")) {
        return(function(object) {
            RET <- pqglht(object)
            m <- coef(object, rhs = TRUE)
            tstat <- switch(object$alternative, 
                            "less" = RET$tstat,
                            "greater" = -RET$tstat,
                            "two.sided" = -abs(RET$tstat))
            C <- object$linfct
            Corder <- C[order(tstat), , drop = FALSE]
            Cm <- m[order(tstat)]
            ms <- maxsets(Corder)
            error <- 0
            p <- sapply(ms, function(x) {
               max(sapply(x, function(s) {
                   object$linfct <- Corder[s, , drop = FALSE]
                   object$rhs <- Cm[s]
                   tmp <- pqglht(object)$pfunction(ifelse(type == "Westfall", 
                                                   "adjusted", "bonferroni"), 
                                                   ...)
                   tmperr <- attr(tmp, "error")
                   if (!is.null(tmperr) && tmperr > error)
                       error <<- tmperr
                   min(tmp)
               }))
            })
            for (i in 2:length(p))
                p[i] <- max(p[i-1], p[i])
            ### <FIXME> what happens in case of ties??? </FIXME> ###
            RET$pvalues <- p[rank(tstat)]
            attr(RET$pvalues, "error") <- error
            RET$type <- type
            class(RET) <- "mtest"
            RET
        })
    }

    ### compute adjustment via p.adjust
    return(function(object) {
        RET <- pqglht(object)
        RET$pvalues <- RET$pfunction(type)
        RET$type <- type
        class(RET) <- "mtest"
        RET
    })
}

adjusted_calpha <- function(...) {
    ret <- function(object, level) {
        pqglht(object)$qfunction(level, adjusted = TRUE, ...)
    }
    attr(ret, "type") <- "adjusted"
    ret
}

univariate_calpha <- function(...) {
    ret <- function(object, level) {
        pqglht(object)$qfunction(level, adjusted = FALSE, ...)
    }
    attr(ret, "type") <- "univariate"
    ret
}

