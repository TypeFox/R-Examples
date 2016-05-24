addbin <- function(y, x, start = NULL, control = list(), allref)
{
    control <- do.call("addreg.control", control)
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    
    nvars <- ncol(x)
    nobs <- NROW(y)
    
    n <- weights <- rep(1, nobs)
    
	bin.id.link <- make.link("identity")
    fam <- binomial(link = bin.id.link)
    eval(fam$initialize)
    
    y1 <- round(n*y)
    y2 <- round(n*(1-y))
    
    Y <- c(y1, y2)
    N <- c(n, n)
    j <- factor(c(rep(1, nobs), rep(2, nobs)))
    
    if (!is.null(start)) {
        start.int <- start[1]
        start.other <- start[-1]
    } else start.int <- start.other <- NULL
    
    if(nvars == 1) {
        X <- data.frame(j)
        names(X) <- "j"
        mono <- FALSE
    } else {
        x.noint <- x[, -1, drop = FALSE]
        x.min <- apply(x.noint, 2, min)
        x.max <- apply(x.noint, 2, max)
        x.res <- sweep(sweep(2*x.noint,MARGIN=2,(x.min + x.max),FUN="-"),MARGIN=2,(x.max - x.min),FUN="/")
        if (!is.null(start)) {
            start.int <- start.int + sum(start.other * (x.min + x.max) / 2)
            start.other <- start.other * (x.max - x.min) / 2
        }
        X <- data.frame(j, rbind(x.res, -x.res), row.names = NULL)
        xnames.temp <- paste("x", seq_len(nvars-1L), sep = "")
        names(X) <- c("j", xnames.temp)
        mono <- FALSE
        termlist <- attr(allref$terms, "term.labels")
        for(term in termlist) {
            if(attr(allref$allref[[term]],"type") == 1) mono <- c(mono, allref$monotonic[term])
            else {
                nlev <- nlevels(factor(allref$data[[term]])) - 1
                mono <- c(mono, rep(TRUE, nlev))
            }
        }
    }

    data.new <- data.frame(Y, N, X)
    if (!is.null(start))
        start.new <- c(start.int, 1 - 2*start.int, start.other)
    else start.new <- NULL

    formula.addpois <- as.formula(paste("Y ~",paste(names(X),collapse=" + ")))
    model.addpois <- addreg(formula.addpois, mono = unname(mono), family = poisson, 
        data = data.new, standard = N, start = start.new, control = control, warn = FALSE, fit = TRUE)
    
    if(nvars == 1) {
        coefs <- model.addpois$coefficients[1]
    } else {
        coefs.int <- model.addpois$coefficients[1]
        coefs.other <- model.addpois$coefficients[-(1:2)]
        coefs <- c(coefs.int - sum(coefs.other*(x.min+x.max)/(x.max-x.min)),
                    2*coefs.other/(x.max-x.min))
    }
    names(coefs) <- xnames
    
	loglik.adj <- sum(lgamma(n+1)) + sum(n*(1-log(n)))
    loglik.bin <- model.addpois$loglik + loglik.adj
	aic.bin <- model.addpois$aic - 2*(1+loglik.adj)
	aic.c.bin <- aic.bin + 2*nvars*(nvars+1)/(nobs - nvars - 1)
	
	wtdmu <- sum(y1) / sum(n)
	nulldev <- sum(fam$dev.resids(y1/n, wtdmu, n))
	nulldf <- nobs - 1
	resdf <- nobs - nvars
	
	list(coefficients = coefs, residuals = model.addpois$residuals[1:nobs],
			fitted.values = model.addpois$fitted.values[1:nobs] / n, rank = model.addpois$rank - 1,
			family = fam, linear.predictors = model.addpois$linear.predictors[1:nobs],
			deviance = model.addpois$deviance, aic = aic.bin, aic.c = aic.c.bin, 
			null.deviance = nulldev, iter = model.addpois$iter, weights = weights,
			prior.weights = n, df.residual = resdf, df.null = nulldf,
			y = y, converged = model.addpois$converged, boundary = model.addpois$boundary,
			loglik = loglik.bin, model.addpois = model.addpois)
}
