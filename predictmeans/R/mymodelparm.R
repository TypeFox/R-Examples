mymodelparm <- function(model, coef., vcov., df, ...) 
UseMethod("mymodelparm")

mymodelparm.default <- function(model, coef. = coef, vcov. = vcov, df = NULL, ...) 
{

    ### extract coefficients and their covariance matrix
    beta <- try(coef.(model))
    if (inherits(beta, "try-error"))
        stop("no ", sQuote("coef"), " method for ",
             sQuote("model"), " found!")

    sigma <- try(vcov.(model))
    if (inherits(sigma, "try-error"))
        stop("no ", sQuote("vcov"), " method for ",
             sQuote("model"), " found!")       
    sigma <- as.matrix(sigma)

    if (any(length(beta) != dim(sigma))) 
        beta = na.omit(beta)
        # stop("dimensions of coefficients and covariance matrix don't match")

    ### determine degrees of freedom
    if (is.null(df)) {
        df <- 0
        ### check if a linear model was supplied
        if (class(model)[1] %in% c("aov", "lm", "glm")) {
            class(model) <- "lm"
            df <- summary(model)$df[2]
        }
		if (class(model)[1] == "gls") {
			dd <- model$dims
			df <- dd[["N"]] - dd[["p"]]
		}
        if (inherits(model, "parm")) df <- model$df
    } else {
        if (df < 0) stop(sQuote("df"), " is not positive")
    }

    ### try to identify non-estimable coefficients
    ### coef.aov removes NAs, thus touch coefficients 
    ### directly
    ocoef <- coef.(model)
    if (inherits(model, "aov")) ocoef <- model$coefficients
    estimable <- rep(TRUE, length(ocoef))
    if (any(is.na(ocoef))) {
        estimable[is.na(ocoef)] <- FALSE
        beta <- ocoef[estimable]
    }

    ### just in case...
    if (length(beta) != ncol(sigma) || nrow(sigma) != sum(estimable))
        stop("could not extract coefficients and covariance matrix from ", 
             sQuote("model"))

    RET <- list(coef = beta, vcov = sigma, df = df, estimable = estimable)
    class(RET) <- "mymodelparm"
    RET
}

mymodelparm.aovlist <- function(model, coef. = coef, vcov. = vcov, df = NULL, ...)
    stop("This function does not support objects of class ", sQuote("aovlist"))

mymodelparm.lme <- function(model, coef. = nlme::fixef, vcov. = vcov, df = NULL, ...)
    mymodelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

mymodelparm.lmerMod <- function(model, coef. = lme4::fixef, vcov. = vcov, df = NULL, ...)
   mymodelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

mymodelparm.glmerMod <- function(model, coef. = lme4::fixef, vcov. = vcov, df = NULL, ...)
   mymodelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)
	
	
	
	