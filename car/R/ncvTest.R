#-------------------------------------------------------------------------------
# Revision history:
# 2009-09-28 by J. Fox (renamed)
# 2012-07-01 Rewritten by S. Weisberg.  The 'data' argument is now gone
# 2013-07-09 Works correctly if data arg is not set in the model
#            works correctly if the formula in 'lm' is an argument
# 2014-11-06 Fixed conflicts with objects in base package. J. Fox
#-------------------------------------------------------------------------------

# score test of nonconstant variance (J. Fox)

ncvTest <- function(model, ...){
    UseMethod("ncvTest")
}

ncvTest.lm <- function(model, var.formula, ...) {
    data <- getCall(model)$data
    model <- if (!is.null(data)){
        data <- eval(data, envir=environment(formula(model)))
        update(model, formula(model), na.action="na.exclude", data=data)
    }
    else update(model, formula(model), na.action="na.exclude")
    sumry <- summary(model)
    residuals <- residuals(model, type="pearson")
    S.sq <- df.residual(model)*(sumry$sigma)^2/sum(!is.na(residuals))
    .U <- (residuals^2)/S.sq
    if (missing(var.formula)) {
        mod <- lm(.U ~ fitted.values(model))
        varnames <- "fitted.values"
        var.formula <- ~ fitted.values
        df <- 1
    }
    else {
        form <- as.formula(paste(".U ~ ", as.character(var.formula)[[2]], sep=""))
        mod <- if(!is.null(data)){
            data$.U <- .U
            lm(form, data=data)
        }  
        else lm(form)
        df <- sum(!is.na(coefficients(mod))) - 1
    }	    
    SS <- anova(mod)$"Sum Sq"
    RegSS <- sum(SS) - SS[length(SS)]
    Chisq <- RegSS/2
    result <- list(formula=var.formula, formula.name="Variance", ChiSquare=Chisq, Df=df, 
        p=pchisq(Chisq, df, lower.tail=FALSE), test="Non-constant Variance Score Test")
    class(result) <- "chisqTest"
    result
}  

ncvTest.glm <- function(model, ...){
    stop("requires lm object")
}

print.chisqTest <- function(x, ...){
    title <- if (!is.null(x$test)) x$test else "Chisquare Test"
    cat(title,"\n")
    if (!is.null(x$formula)) cat(x$formula.name, 
        "formula:", as.character(x$formula), "\n")
    cat("Chisquare =", x$ChiSquare,"   Df =", x$Df,
        "    p =", x$p, "\n")
    invisible(x)
}
