#' @export

model.matrix.rdd_data <- function(object, covariates = NULL, order = 1, bw = NULL, slope = c("separate", "same"), covar.opt = list(strategy = c("include", 
    "residual"), slope = c("same", "separate"), bw = NULL), covar.strat = c("include", "residual"), ...) {
    
    checkIsRDD(object)
    rdd_object <- object
    type <- getType(object)
    
    if (!missing(covar.strat)) 
        warning("covar.strat is (soon) deprecated arg!")
    
    slope <- match.arg(slope)
    covar.strat <- match.arg(covar.opt$strategy, choices = c("include", "residual"))
    covar.slope <- match.arg(covar.opt$slope, choices = c("same", "separate"))
    
    cutpoint <- getCutpoint(rdd_object)
    if (!is.null(covariates) & !hasCovar(rdd_object)) 
        stop("Arg 'covariates' was specified, but no covariates found in 'rdd_object'.")
    
    ## Construct data
    dat <- as.data.frame(rdd_object)
    
    dat_step1 <- dat[, c("y", "x")]
    dat_step1$x <- dat_step1$x - cutpoint
    
    L <- ifelse(dat_step1$x >= 0, 1, 0)
    dat_step1$D <- if (type == "Sharp") 
        L else object$z
    
    if (order > 0) {
        polys <- poly(dat_step1$x, degree = order, raw = TRUE)
        colnames(polys) <- paste("x", 1:order, sep = "^")
        dat_step1 <- cbind(dat_step1[, c("y", "D")], polys)
        if (slope == "separate") {
            polys2 <- polys * L
            colnames(polys2) <- paste(colnames(polys), "right", sep = "_")
            dat_step1 <- cbind(dat_step1, polys2)
        }
    } else {
        dat_step1$x <- NULL
    }
    
    ## Covariates
    if (!is.null(covariates)) {
        covar <- getCovar(rdd_object)
        formu.cova <- covariates
        
        if (grepl("\\.", formu.cova)) 
            formu.cova <- paste(colnames(covar), collapse = " + ")
        if (covar.slope == "separate") {
            formu.cova <- paste(formu.cova, "+", paste("D*(", formu.cova, ")", sep = ""), sep = " ")
            covar$D <- dat_step1$D
        }
        
        formula.cova <- as.formula(paste("~", formu.cova))
        mf <- model.frame(formula.cova, covar, na.action = na.pass)
        M_covar <- model.matrix(formula.cova, data = mf)
        
        if (covar.strat == "residual") {
            M_covar <- data.frame(y = dat_step1$y, M_covar)
            first_stage <- lm(y ~ ., data = M_covar)  ## regress y on covariates only
            dat_step1$y <- residuals(first_stage)  ## change in original data
        } else {
            rem <- switch(covar.slope, separate = "^D$|(Intercept)", same = "(Intercept)")
            M_covar <- M_covar[, -grep(rem, colnames(M_covar)), drop = FALSE]
            dat_step1 <- cbind(dat_step1, M_covar)  ## add covar as regressors
        }
    }
    
    ## Colnames cleaning
    colnames(dat_step1) <- gsub("x\\^1", "x", colnames(dat_step1))
    
    ## 
    if (type == "Fuzzy") 
        dat_step1$ins <- L
    
    ## return results:
    dat_step1
} 
