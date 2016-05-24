feature.test <- function(x, y, B=100,
                         type.measure="deviance", s="lambda.min",
                         keeplambda=FALSE,
                         olsestimates=TRUE,
                         penalty.factor = rep(1, nvars),
                         alpha=1,
                         control=list(trace=FALSE, maxcores=24), ...) {

#    require(parallel)
#    require(glmnet)

    # Suppress warnings while running
    warn <- options()$warn
    options(warn=-1)
    
    
    # Parse the control options
    con <- list(trace=FALSE)
    con[names(control)] <- control
    control <- con

    sfixed <- ifelse(is.numeric(s), TRUE, FALSE)
    lambda <- NA

    nobs <- nrow(x)
    nvars <- ncol(x)

    # Starts by scaling both the X and the Y variable
    # This should save us a wee bit of speed later on
    x <- scale(x)
    y <- y / sd(y)


    # Do we need to estimate lambda?
    if (!sfixed) {
        o <- cv.glmnet(x, y, standardize=FALSE, type.measure=type.measure, penalty.factor=penalty.factor, ..., grouped=FALSE, pmax = min(nobs-1, nvars))
        lambda <- as.numeric(o[s])
    } else {
        lambda <- s
        keeplambda <- TRUE       
    }
    
    # We assume that the intercept is automatically in the model
    o <- glmnet(x, y, standardize=FALSE, penalty.factor=penalty.factor, ..., pmax = min(nobs-1, nvars))
    # Compute the deviance explained (using linear interpolation)
#    o.devexp <- approx(o$lambda, o$dev.ratio, lambda)$y

    # Now find the non-zero predictors
    o.coef <- coef(o, s=lambda)[-1]    # We remove the intercept as that is always the first predictor
    o.select <- which(o.coef != 0)     # Contains the index of variables found to be non-zero (apart from intercept)
    o.non.zero <- o.coef[o.select]     # Value of coefficients actually non-zero
    o.nselected  <- length(o.non.zero) # Number of non-zero predictors
    o.sigmalasso <- sqrt(sum( (y - predict(o, newx=x, s=lambda))^2)/(nobs - o.nselected -1))   # Estimate the residual variance
    o.coef.lasso <- o.non.zero
    names(o.coef.lasso) <- o.select
#    o.lasso.order <- order(abs(o.non.zero), decreasing=TRUE)
#    o.lasso.coef <- abs(o.coef.lasso[o.lasso.order])

    
    # Compute ols estimates (if we found something)
    if (o.nselected>0) {
        o.lm <- lm(y ~ x[,o.select])
        o.summarylm <- summary(o.lm)
#        print(o.summarylm)
#        o.rsqr <- o.summarylm$r.squared
        o.fullp <- -pf(o.summarylm$fstatistic[1], o.summarylm$fstatistic[2], o.summarylm$fstatistic[3], lower.tail=FALSE, log.p=TRUE)
        o.tstat <- coef(o.summarylm)[-1,3]
        o.betacoef <- coef(o.summarylm)[-1,1]
        names(o.tstat) <- o.select
        o.sigmaols <- o.summarylm$sigma

#        o.ols.order <- order(abs(o.tstat), decreasing=TRUE)
#        o.ols.coef <- abs(o.tstat)[o.ols.order]    # Ordered t-test statistics from ols
    }
    else {
        o.betacoef <- 0
        o.sigmalasso <- 0
        o.sigmaols <- 0
        o.coef.lasso <- 0
#        o.ols.coef <- NA
        o.ols.order <- NA
        o.fullp <- 0
        o.tstat <- 0
    }


#    print(o.fullp)
    
#    cat("Lasso results\n")
#    cat("  Selected:\n ")
#    print(o.select)
#    cat("  Beta-hat according to size")
#    print(o.lasso.coef)
#    print(o.lasso.order)

#    cat("OLS\n")
#    cat("  Beta-hat according to size")
#    print(o.ols.coef)
#    print(o.ols.order)

    scaling.factor <- o.coef.lasso * o.sigmaols / (abs(o.betacoef) * o.sigmalasso )
    o.lasso.teststat <- scaling.factor*o.tstat

    o.lasso.torder <- order(abs(o.lasso.teststat), decreasing=TRUE)


    o.lasso.max <- max(c(0, abs(o.lasso.teststat)))
    o.ols.max <- max(c(0, abs(o.tstat)))

    # Run the simulations under the NULL
    if (o.nselected > 0) {
        simnull <- mclapply(1:B, function(iii) {

            py <- sample(y)
            cat(".")

            if (!keeplambda) {
                to <- cv.glmnet(x, py, standardize=FALSE, type.measure=type.measure, penalty.factor=penalty.factor, ..., pmax = min(nobs-1, nvars), grouped=FALSE)
                lambda <- as.numeric(to[s])
            }

            # We assume that the intercept is automatically in the model
            to <- glmnet(x, py, standardize=FALSE, penalty.factor=penalty.factor, ..., pmax = min(nobs-1, nvars))
            # Compute the deviance explained (using linear interpolation)
#            to.devexp <- approx(to$lambda, to$dev.ratio, lambda)$y

            # Now find the non-zero predictors
            to.coef <- coef(to, s=lambda)[-1]    # We remove the intercept as that is always the first predictor
            to.select <- which(to.coef != 0)     # Contains the index of variables found to be non-zero (apart from intercept)
            to.non.zero <- to.coef[to.select]     # Value of coefficients actually non-zero
            to.nselected  <- length(to.non.zero) # Number of non-zero predictors
            to.sigmalasso <- sqrt(sum( (py - predict(to, newx=x, s=lambda))^2)/(nobs - to.nselected -1))   # Estimate the residual variance
            to.coef.lasso <- to.non.zero
            names(to.coef.lasso) <- to.select
#            to.lasso.order <- order(abs(to.non.zero), decreasing=TRUE)
#            to.lasso.coef <- abs(to.coef.lasso[to.lasso.order])

            # Compute ols estimates (but only if we found something)
            if (to.nselected>0) {
                to.lm <- lm(py ~ x[,to.select])
                to.summarylm <- summary(to.lm)
#                print(to.summarylm)
#                to.rsqr <- to.summarylm$r.squared
                to.fullp <- -pf(to.summarylm$fstatistic[1], to.summarylm$fstatistic[2], to.summarylm$fstatistic[3], lower.tail=FALSE, log.p=TRUE)
                to.tstat <- coef(to.summarylm)[-1,3]
                to.betacoef <- coef(to.summarylm)[-1,1]
                names(to.tstat) <- to.select
                to.sigmaols <- to.summarylm$sigma

#                to.ols.order <- order(abs(to.tstat), decreasing=TRUE)
#                to.ols.coef <- abs(to.tstat)[to.ols.order]    # Ordered t-test statistics from ols
            } else {
#                to.ols.order <- NA
#                to.ols.coef <- NA
                to.sigmaols <- 0
                to.betacoef <- 0
                to.tstat <- 0
                to.fullp <- 0
            }

#            cat("Lasso results\n")
#            cat("  Selected:\n ")
#            print(to.select)
#            cat("  Beta-hat according to size")
#            print(to.lasso.coef)
#            print(to.lasso.order)
#            cat("OLS\n")
#            cat("  Beta-hat according to size")
#            print(to.ols.coef)
#            print(to.ols.order)
#            cat("-------\n")

#            print(c(to.coef.lasso, to.sigmaols, (abs(to.betacoef)), to.sigmalasso ))
#            print(to.tstat)
            scaling.factor <- to.coef.lasso * to.sigmaols / (abs(to.betacoef) * to.sigmalasso )
            to.lasso.teststat <- scaling.factor*max(abs(to.tstat))

            to.lasso.max <- max(c(0, to.lasso.teststat))
            to.ols.max <- max(c(0, abs(to.tstat)))
#            cat("====\n")

#            print(c(to.lasso.max, to.ols.max, to.fullp))
            c(to.lasso.max, to.ols.max, to.fullp)
            
        }, mc.cores=min(detectCores(), control$maxcores)  )

        # Make the summary statistics and evaluate
        dist.lasso <- sapply(simnull, function(i) {i[[1]]})
        dist.ols   <- sapply(simnull, function(i) {i[[2]]})        
        dist.fullp <- sapply(simnull, function(i) {i[[3]]})

        p.maxlasso <- sum(dist.lasso >= o.lasso.max)/B
        p.maxols   <- sum(dist.ols >= o.ols.max)/B
        p.full     <- sum(dist.fullp >= o.fullp)/B

        
        
    } else {
        # None selected from the original data

        p.maxlasso <- NA
        p.maxols   <- NA
        p.full     <- NA

        
    }
    

    options(warn=warn)
    
    list(
        p.full=p.full,
        ols.selected=o.select[o.ols.order],
        p.maxols=p.maxols,
        lasso.selected=o.select[o.lasso.torder],            
        p.maxlasso=p.maxlasso,
        lambda.orig=lambda,
        B=B)
    
}


