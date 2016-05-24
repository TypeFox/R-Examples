# this is a test summary for using the vcov matrix
# it should replace  he old summary.gamlss
#---------------------------------------------------------------------------------------
summary.gamlss<- function (object,
                           type = c("vcov", "qr"), robust = FALSE,
                           save = FALSE, 
                    hessian.fun = c("R", "PB"), 
                         digits = max(3, getOption("digits") - 3),
                           ...) 
{
type <- match.arg(type)
## if use the vcov for se's
if (type=="vcov")
 {
## this will be moved  
      covmat <- try(suppressWarnings(vcov(object, type="all", robust=robust,  hessian.fun = hessian.fun)), silent = TRUE) 
             if (any(class(covmat)%in%"try-error"||any(is.na(covmat$se))))
                 { 
                   warning(paste("summary: vcov has failed, option qr is used instead\n"))
                   type <- "qr"
                 }
## if try fails call old summary
 }      
ifWarning <- rep(FALSE, length(object$parameters))# to create warnings
if (type=="vcov")#  type vcov --------------------------------------------------
 {
        coef <- covmat$coef
          se <- covmat$se
      tvalue <- coef/se
      pvalue <-  2 * pt(-abs(tvalue), object$df.res)  #if (est.disp) 2 * pt(-abs(tvalue), df.r) else   2 * pnorm(-abs(tvalue)) 
  coef.table <- cbind(coef, se, tvalue, pvalue)
  dimnames(coef.table) <- list(names(coef), c("Estimate" , "Std. Error" ,"t value","Pr(>|t|)"))  
#now the table contains all informatio I need 
#browser()
#printCoefmat(coef.table, digits = digits, signif.stars = TRUE)
# to print use   coef.table[1:3,] 
      cat("*******************************************************************")
    cat("\nFamily: ", deparse(object$family), "\n") 
    cat("\nCall: ", deparse(object$call),  "\n", fill=TRUE)
    cat("Fitting method:", deparse(object$method), "\n\n") 
    est.disp <- FALSE
       # df.r <- object$noObs - object$mu.df
#================ mu ESTIMATES ========================
    if ("mu"%in%object$parameters)   
    {
      
      ifWarning[1]  <- (!is.null(unlist(attr(terms(formula(object, "mu"), 
                  specials = .gamlss.sm.list), "specials")))) 
        if (object$mu.df != 0)   
        {          
                        pm <- object$mu.qr$rank 
                        p1 <- 1:pm
            cat("-------------------------------------------------------------------\n")
            cat("Mu link function: ", object$mu.link)
            cat("\n")
            cat("Mu Coefficients:")
            if (is.character(co <- object$contrasts)) 
                cat("  [contrasts: ", apply(cbind(names(co), co), 1, 
                    paste, collapse = "="), "]")
            cat("\n")
            printCoefmat(coef.table[p1,,drop=FALSE], digits = digits, signif.stars = TRUE)
     #       print.default(coef.table[p1,], digits = digits, print.gap = 2, quote = FALSE)
            cat("\n")
        }
        else
            if(object$mu.fix == TRUE) 
            {
            #warning(paste("mu parameter is fixed"))
            cat("-------------------------------------------------------------------\n")
            cat("Mu parameter is fixed \n")
            if (all(object$mu.fv == object$mu.fv[1]))
                cat("Mu = ", object$mu.fv[1], "\n")
            else
                cat("Mu is equal with the vector (", object$mu.fv[1], ",",object$mu.fv[2], ",",object$mu.fv[3], ",",object$mu.fv[4], ", ...) \n")
            }
    }
    if ("sigma"%in%object$parameters) 
    {
      ifWarning[2]  <- (!is.null(unlist(attr(terms(formula(object, "sigma"), 
                     specials = .gamlss.sm.list), "specials")))) 
         if (object$sigma.df != 0)   
        {
                     ps <- object$sigma.qr$rank 
                     p1 <- (pm+1):(pm+ps)
        #cat("\n")
        cat("-------------------------------------------------------------------\n")
        cat("Sigma link function: ", object$sigma.link)
        cat("\n")
        cat("Sigma Coefficients:")
        cat("\n")    
        printCoefmat(coef.table[p1,, drop=FALSE], digits = digits, signif.stars = TRUE)
       # print.default(coef.table[p1,], digits = digits, print.gap = 2,  quote = FALSE)
        cat("\n")
        }
        else 
            if(object$sigma.fix == TRUE) 
            {
            cat("-------------------------------------------------------------------\n")
            cat("Sigma parameter is fixed")
            cat("\n")
            if (all(object$sigma.fv == object$sigma.fv[1]))
                cat("Sigma = ", object$sigma.fv[1], "\n")
            else
                cat("Sigma is equal with the vector (", object$sigma.fv[1], ",",object$sigma.fv[2], ",",object$sigma.fv[3], ",",object$sigma.fv[4], ", ...) \n")
            }
    }
    if ("nu"%in%object$parameters) 
    {
      ifWarning[3]  <- (!is.null(unlist(attr(terms(formula(object, "nu"), 
                      specials = .gamlss.sm.list), "specials")))) 
         if (object$nu.df != 0)   
        {
                   
                     pn <- object$nu.qr$rank 
                     p1 <- (pm+ps+1):(pm+ps+pn)
        cat("-------------------------------------------------------------------\n")
        cat("Nu link function: ", object$nu.link,"\n")
        cat("Nu Coefficients:")
        cat("\n")
        printCoefmat(coef.table[p1,, drop=FALSE], digits = digits, signif.stars = TRUE)
       # print.default(coef.table[p1,], digits = digits, print.gap = 2,  quote = FALSE)
        cat("\n")
        }
        else
            if(object$nu.fix == TRUE) {
            cat("-------------------------------------------------------------------\n")
            cat("Nu parameter is fixed")
            cat("\n")
            if (all(object$nu.fv == object$nu.fv[1]))
                cat("Nu = ", object$nu.fv[1], "\n")
            else
                cat("Nu is equal with the vector (", object$nu.fv[1], ",",object$nu.fv[2], ",",object$nu.fv[3], ",",object$nu.fv[4], ", ...) \n")
            }
    }
    if ("tau"%in%object$parameters) 
    {
      ifWarning[4]  <- (!is.null(unlist(attr(terms(formula(object, "tau"), 
                     specials = .gamlss.sm.list), "specials")))) 
         if (object$tau.df != 0)   
        {
                    
                     pt <- object$tau.qr$rank 
                     p1 <- (pm+ps+pn+1):(pm+ps+pn+pt)
        cat("-------------------------------------------------------------------\n")
        cat("Tau link function: ", object$tau.link,"\n")
        cat("Tau Coefficients:")
        cat("\n")
        printCoefmat(coef.table[p1,, drop=FALSE], digits = digits, signif.stars = TRUE)
       # print.default(coef.table[p1,], digits = digits, print.gap = 2,  quote = FALSE)
        cat("\n")
        }
        else
            if(object$tau.fix == TRUE) {
            cat("-------------------------------------------------------------------\n")
            cat("Tau parameter is fixed")
            cat("\n")
            if (all(object$tau.fv == object$tau.fv[1]))
                cat("Tau = ", object$tau.fv[1], "\n")
            else
                cat("Tau is equal with the vector (", object$tau.fv[1], ",",object$tau.fv[2], ",",object$tau.fv[3], ",",object$tau.fv[4], ", ...) \n")
            }
    }
if (any(ifWarning))
{
  cat("-------------------------------------------------------------------\n")
  cat("NOTE: Additive smoothing terms exist in the formulas: \n")
  cat(" i) Std. Error for smoothers are for the linear effect only. \n")
  cat("ii) Std. Error for the linear terms maybe are not accurate. \n")
}
   cat("-------------------------------------------------------------------\n")
   cat("No. of observations in the fit: ", object$noObs, "\n")
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
if (type=="qr")#     TYPE qr ---------------------------------------------------
  {
# local function definition
#-------------------------------------------------------------------------------
        estimatesgamlss<-function (object,Qr, p1, coef.p, est.disp , df.r, digits = max(3, getOption("digits") - 3),
                                               covmat.unscaled , ...)
          {
                             #covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
             dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
                                covmat <- covmat.unscaled #in glm is=dispersion * covmat.unscaled, but here is already multiplied by the dispersion
                                var.cf <- diag(covmat)
                                 s.err <- sqrt(var.cf)
                                tvalue <- coef.p/s.err
                                    dn <- c("Estimate", "Std. Error")
                          if (!est.disp) 
                             {
                                            pvalue <- 2 * pnorm(-abs(tvalue))
                                        coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
                              dimnames(coef.table) <- list(names(coef.p), c(dn, "z value",
                                                    "Pr(>|z|)"))
                             }
                           else if (df.r > 0) 
                             {
                                           pvalue <- 2 * pt(-abs(tvalue), df.r)
                                       coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
                             dimnames(coef.table) <- list(names(coef.p), c(dn, "t value","Pr(>|t|)"))
                             }
                           else 
                           {
                                          coef.table <- cbind(coef.p, Inf)
                                dimnames(coef.table) <- list(names(coef.p), dn)
                           }
                          return(coef.table)
          }
##------------------------------------------------------------------------------
# here for the proper function
##------------------------------------------------------------------------------
      dispersion <- NULL
 #     digits <- max(3, getOption("digits") - 3)
      cat("*******************************************************************")
    cat("\nFamily: ", deparse(object$family), "\n") 
      cat("\nCall: ", deparse(object$call),  "\n", fill=TRUE)
    cat("Fitting method:", deparse(object$method), "\n\n") 
    est.disp <- FALSE
        df.r <- object$noObs - object$mu.df
#================ mu ESTIMATES ========================
    if ("mu"%in%object$parameters)   
    {
      ifWarning[1]  <- (!is.null(unlist(attr(terms(formula(object, "mu"), 
                             specials = .gamlss.sm.list), "specials")))) 
        if (object$mu.df != 0)   
        {
            Qr <- object$mu.qr 
            df.r <- object$noObs - object$mu.df
            ## this should be taken out sinse there is no dispersion MS Thursday, August 17, 2006
            if (is.null(dispersion)) 
                dispersion <- if (any(object$family == c("PO", 
                                    "BI", "EX", "P1"))) 
                                    1
                            else if (df.r > 0) {
                                                est.disp <- TRUE
                                                if (any(object$weights == 0)) 
                                                        warning(paste("observations with zero weight", 
                                                        "not used for calculating dispersion"))
                                                }
            else Inf   
            #---------------------------------------------------end of taken out 
                        p <- object$mu.df #  object$rank
                        p1 <- 1:(p-object$mu.nl.df)
            covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
            mu.coef.table <- estimatesgamlss(object=object,Qr=object$mu.qr, p1=p1, 
                                        coef.p=object$mu.coefficients[Qr$pivot[p1]], 
                                        est.disp =est.disp, df.r=df.r,
                                        covmat.unscaled =covmat.unscaled )
            cat("-------------------------------------------------------------------\n")
            cat("Mu link function: ", object$mu.link)
            cat("\n")
            cat("Mu Coefficients:")
            if (is.character(co <- object$contrasts)) 
                cat("  [contrasts: ", apply(cbind(names(co), co), 1, 
                    paste, collapse = "="), "]")
            cat("\n")
            printCoefmat(mu.coef.table, digits = digits, signif.stars = TRUE)
           # print.default(mu.coef.table, digits = digits, print.gap = 2, quote = FALSE)
            cat("\n")
        }
        else
            if(object$mu.fix == TRUE) 
            {
            #warning(paste("mu parameter is fixed"))
            cat("-------------------------------------------------------------------\n")
            cat("Mu parameter is fixed")
            cat("\n")
            if (all(object$mu.fv == object$mu.fv[1]))
                cat("Mu = ", object$mu.fv[1], "\n")
            else
                cat("Mu is equal with the vector (", object$mu.fv[1], ",",object$mu.fv[2], ",",object$mu.fv[3], ",",object$mu.fv[4], ", ...) \n")
            }
        coef.table <- mu.coef.table
    }
    else
    {
        if (df.r > 0) {
                        est.disp <- TRUE
                        if (any(object$weights == 0)) 
                        warning(paste("observations with zero weight", 
                        "not used for calculating dispersion"))
                        }
    }
    if ("sigma"%in%object$parameters) 
    {
      ifWarning[2]  <- (!is.null(unlist(attr(terms(formula(object, "sigma"), 
                        specials = .gamlss.sm.list), "specials")))) 
         if (object$sigma.df != 0)   
        {
                      Qr <- object$sigma.qr 
                    df.r <- object$noObs - object$sigma.df
                       p <- object$sigma.df 
                      p1 <- 1:(p-object$sigma.nl.df)
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
        sigma.coef.table<-estimatesgamlss(object=object,Qr=object$sigma.qr, p1=p1, 
                                       coef.p=object$sigma.coefficients[Qr$pivot[p1]], 
                                       est.disp =est.disp, df.r=df.r,
                                       covmat.unscaled =covmat.unscaled )
        #cat("\n")
        cat("-------------------------------------------------------------------\n")
        cat("Sigma link function: ", object$sigma.link)
        cat("\n")
        cat("Sigma Coefficients:")
        cat("\n")
        printCoefmat(sigma.coef.table, digits = digits, signif.stars = TRUE)
        #print.default(sigma.coef.table, digits = digits, print.gap = 2,  quote = FALSE)
        cat("\n")
        }
        else 
            if(object$sigma.fix == TRUE) {
            cat("-------------------------------------------------------------------\n")
            cat("Sigma parameter is fixed")
            cat("\n")
            if (all(object$sigma.fv == object$sigma.fv[1]))
                cat("Sigma = ", object$sigma.fv[1], "\n")
            else
                cat("Sigma is equal with the vector (", object$sigma.fv[1], ",",object$sigma.fv[2], ",",object$sigma.fv[3], ",",object$sigma.fv[4], ", ...) \n")
            }
         coef.table <- rbind(mu.coef.table, sigma.coef.table)
    }
    if ("nu"%in%object$parameters) 
    {
      ifWarning[3]  <- (!is.null(unlist(attr(terms(formula(object, "nu"), 
                             specials = .gamlss.sm.list), "specials")))) 
         if (object$nu.df != 0)   
        {
                     Qr <- object$nu.qr 
                   df.r <- object$noObs - object$nu.df
                      p <- object$nu.df 
                     p1 <- 1:(p-object$nu.nl.df)
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
          nu.coef.table <- estimatesgamlss(object=object,Qr=object$nu.qr, p1=p1, 
                                       coef.p=object$nu.coefficients[Qr$pivot[p1]], 
                                       est.disp =est.disp, df.r=df.r,
                                       covmat.unscaled =covmat.unscaled )
        cat("-------------------------------------------------------------------\n")
        cat("Nu link function: ", object$nu.link,"\n")
        cat("Nu Coefficients:")
        cat("\n")
        printCoefmat(nu.coef.table, digits = digits, signif.stars = TRUE)
        #print.default(nu.coef.table, digits = digits, print.gap = 2,  quote = FALSE)
        cat("\n")
        }
        else
            if(object$nu.fix == TRUE) {
            cat("-------------------------------------------------------------------\n")
            cat("Nu parameter is fixed")
            cat("\n")
            if (all(object$nu.fv == object$nu.fv[1]))
                cat("Nu = ", object$nu.fv[1], "\n")
            else
                cat("Nu is equal with the vector (", object$nu.fv[1], ",",object$nu.fv[2], ",",object$nu.fv[3], ",",object$nu.fv[4], ", ...) \n")
            }
         coef.table <- rbind(mu.coef.table, sigma.coef.table, nu.coef.table)     
    }

    if ("tau"%in%object$parameters) 
   {
      ifWarning[4]  <- (!is.null(unlist(attr(terms(formula(object, "tau"), 
                             specials = .gamlss.sm.list), "specials")))) 
         if (object$tau.df != 0)   
        {
                     Qr <- object$tau.qr 
                   df.r <- object$noObs - object$tau.df
                      p <- object$tau.df 
                     p1 <- 1:(p-object$tau.nl.df)
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
         tau.coef.table <- estimatesgamlss(object=object,Qr=object$tau.qr, p1=p1, 
                                       coef.p=object$tau.coefficients[Qr$pivot[p1]], 
                                       est.disp =est.disp, df.r=df.r,
                                       covmat.unscaled =covmat.unscaled )
        cat("-------------------------------------------------------------------\n")
        cat("Tau link function: ", object$tau.link,"\n")
        cat("Tau Coefficients:")
        cat("\n")
        printCoefmat(tau.coef.table, digits = digits, signif.stars = TRUE)
        #print.default(tau.coef.table, digits = digits, print.gap = 2,  quote = FALSE)
        cat("\n")
        }
        else
            if(object$tau.fix == TRUE) {
            cat("-------------------------------------------------------------------\n")
            cat("Tau parameter is fixed")
            cat("\n")
            if (all(object$tau.fv == object$tau.fv[1]))
                cat("Tau = ", object$tau.fv[1], "\n")
            else
                cat("Tau is equal with the vector (", object$tau.fv[1], ",",object$tau.fv[2], ",",object$tau.fv[3], ",",object$tau.fv[4], ", ...) \n")
            }
         coef.table <- rbind(mu.coef.table, sigma.coef.table, nu.coef.table, tau.coef.table)    
   }
if (any(ifWarning))
{
  cat("-------------------------------------------------------------------\n")
  cat("NOTE: Additive smoothing terms exist in the formulas: \n")
  cat(" i) Std. Error for smoothers are for the linear effect only. \n")
  cat("ii) Std. Error for the linear terms may not be reliable. \n")
}

   cat("-------------------------------------------------------------------\n")
   cat("No. of observations in the fit: ", object$noObs, "\n")
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
 if ( save == TRUE)
    {
      out <- as.list(environment())
       return(out)
  }
invisible(coef.table)
}
 
#==============================================================================
