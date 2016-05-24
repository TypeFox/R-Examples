vf <- function (ord, vars, nperm = 100) 
{
    vfcalc <- function(ord, vars) {
        lm.list <- apply(vars, 2, function(x, ord) lm(x ~ ord), 
            ord = ord)
        coef.m <- sapply(lm.list, function(x) unlist(x$coefficients))
        scores <- atan(coef.m[3, ]/coef.m[2, ])
        scores <- t(rbind(cos(scores), sin(scores)))
        scores <- abs(scores)
        cor.m <- t(cor(ord, vars))
        cor.m <- sign(cor.m)
        scores <- scores * cor.m
        coef.m <- lapply(lm.list, summary)
        r <- sapply(coef.m, function(x) sqrt(unlist(x$r.squared)))
        list(scores = scores, r = r)
    }
    ord <- as.matrix(ord)
    if(is.vector(vars)) {
    	vars <- matrix(vars, ncol=1)
	colnames(vars) <- "var1"
    }
    else
    	vars <- as.matrix(vars)
    nvars <- ncol(vars)
    if (any(is.na(vars))) {
        warning("NA values in variables will be removed\n")
        naindex <- apply(vars, 1, function(x) any(is.na(x)))
        ord <- ord[!naindex, ]
        vars <- vars[!naindex, ]
    }
    vf1 <- vfcalc(ord, vars)
    if (nperm > 0) {
        how.many <- rep(nrow(ord), nperm - 1)
        perm.ord <- lapply(how.many, function(x) sample(1:x))
        r.list <- sapply(perm.ord, function(x, ord, vars, f) f(ord[x, 
            ], vars)$r, ord = ord, vars = vars, f = vfcalc)
	if(nvars == 1) 
	    r.list <- matrix(r.list, nrow=1)
        r.list <- cbind(vf1$r, r.list)
        pval <- apply(r.list, 1, function(x, nperm) length(x[x >= 
            x[1]])/nperm, nperm = nperm)
    }
    else pval <- rep(0, ncol(vars))
    vfres <- cbind(vf1$scores, vf1$r, pval)
    dimnames(vfres)[[1]] <- dimnames(vars)[[2]]
    dimnames(vfres)[[2]] <- c(1:ncol(ord), "r", "pval")
    class(vfres) <- "vf"
    vfres
}

