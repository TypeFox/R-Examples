# KA+MS last change Friday, Thursday, January 2, 2003 at 17:08
summary.nlgamlss<- function (object, dispersion = NULL,  ...) 
{

# here for the proper function
      digits <- max(3, getOption("digits") - 3)
      cat("*******************************************************************")
    cat("\nFamily: ", deparse(object$family), "\n") 
    cat("\nCall: ", deparse(object$call),  "\n", fill=TRUE)
    cat("Fitting method:", deparse(object$method), "\n\n") 
    est.disp <- FALSE
        df.r <- object$noObs - object$mu.df
#================ mu ESTIMATES ========================
if ("mu"%in%object$parameters)   
    {
        if (object$mu.df != 0)   
        { 
                 tvalue <- coef(object,"mu")/object$mu.se
                     dn <- c("Estimate", "Std. Error", "t-value", "p-value")
                 pvalue <- 2 * pnorm(-abs(tvalue))
             mu.coef.table <- cbind(coef(object,"mu"),object$mu.se, tvalue, pvalue)
             colnames(mu.coef.table) <- dn
            cat("-------------------------------------------------------------------\n")
            cat("Mu link function: ", object$mu.link)
            cat("\n")
            cat("Mu Coefficients:")
            if (is.character(co <- object$contrasts)) 
                cat("  [contrasts: ", apply(cbind(names(co), co), 1, 
                    paste, collapse = "="), "]")
            cat("\n")
            print.default(mu.coef.table, digits = digits, print.gap = 2, quote = FALSE)
            cat("\n")
        }
        else
        {
            cat("-------------------------------------------------------------------\n")
            cat("Mu parameter is fixed")
            cat("\n")
            
        }
    }
 if ("sigma"%in%object$parameters) 
    {
         if (object$sigma.df != 0)   
        {
                          tvalue <- coef(object,"sigma")/object$sigma.se
                     dn <- c("Estimate", "Std. Error", "t-value", "p-value")
                 pvalue <- 2 * pnorm(-abs(tvalue))
       sigma.coef.table <- cbind(coef(object,"sigma"),object$sigma.se, tvalue, pvalue)
    colnames(sigma.coef.table) <- dn
            cat("-------------------------------------------------------------------\n")
            cat("Sigma link function: ", object$sigma.link)
            cat("\n")
            cat("Migma Coefficients:")
            if (is.character(co <- object$contrasts)) 
                cat("  [contrasts: ", apply(cbind(names(co), co), 1, 
                    paste, collapse = "="), "]")
            cat("\n")
            print.default(sigma.coef.table, digits = digits, print.gap = 2, quote = FALSE)
            cat("\n")
        }
            else
        {
            cat("-------------------------------------------------------------------\n")
            cat("Sigma parameter is fixed")
            cat("\n")
            
        }
    }
    if ("nu"%in%object$parameters) 
    {
        if (object$nu.df != 0)   
        { 
                 tvalue <- coef(object,"nu")/object$nu.se
                     dn <- c("Estimate", "Std. Error", "t-value", "p-value")
                 pvalue <- 2 * pnorm(-abs(tvalue))
             nu.coef.table <- cbind(coef(object,"nu"),object$nu.se, tvalue, pvalue)
             colnames(nu.coef.table) <- dn
            cat("-------------------------------------------------------------------\n")
            cat("Nu link function: ", object$nu.link)
            cat("\n")
            cat("Nu Coefficients:")
            if (is.character(co <- object$contrasts)) 
                cat("  [contrasts: ", apply(cbind(names(co), co), 1, 
                    paste, collapse = "="), "]")
            cat("\n")
            print.default(nu.coef.table, digits = digits, print.gap = 2, quote = FALSE)
            cat("\n")
        }
        else
        {
            cat("-------------------------------------------------------------------\n")
            cat("Nu parameter is fixed")
            cat("\n")
            
        }
    }
    if ("tau"%in%object$parameters) 
   {
        if (object$tau.df != 0)   
        { 
                 tvalue <- coef(object,"tau")/object$tau.se
                     dn <- c("Estimate", "Std. Error", "t-value", "p-value")
                 pvalue <- 2 * pnorm(-abs(tvalue))
             tau.coef.table <- cbind(coef(object,"tau"),object$tau.se, tvalue, pvalue)
             colnames(tau.coef.table) <- dn
            cat("-------------------------------------------------------------------\n")
            cat("Tau link function: ", object$tau.link)
            cat("\n")
            cat("Tau Coefficients:")
            if (is.character(co <- object$contrasts)) 
                cat("  [contrasts: ", apply(cbind(names(co), co), 1, 
                    paste, collapse = "="), "]")
            cat("\n")
            print.default(tau.coef.table, digits = digits, print.gap = 2, quote = FALSE)
            cat("\n")
        }
        else
        {
            cat("-------------------------------------------------------------------\n")
            cat("Tau parameter is fixed")
            cat("\n")
            
        }
    }
   cat("-------------------------------------------------------------------\n")
   cat("No. of observations in the fit: ", object$N, "\n")
   cat("Degrees of Freedom for the fit: ", object$df.fit)
   cat("\n")
   cat("      Residual Deg. of Freedom: ", object$df.residual, "\n")
   cat("                      at cycle: ", object$iter, "\n \n")
   cat("Global Deviance:    ", object$G.deviance,#format(signif(object$G.deviance, digits)), 
        "\n            AIC:    ",object$aic, #format(signif(object$aic, digits)), 
        "\n            SBC:    ",object$sbc, "\n") #format(signif(object$sbc, digits)), "\n")
    cat("*******************************************************************")
    cat("\n")
}

#==============================================================================
vcov.nlgamlss <- function (object, ...) 
{ 
    object$cov
}
#==============================================================================
