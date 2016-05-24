## last modified 2011-10-18

anova.mix <- function (object, mixobj2, ...) 
{
    mixobj1 <- object
    mixanova <- NULL
    if (inherits(mixobj1, "mix")) {
        hypo <- FALSE
        if (!missing(mixobj2)) {
            if (inherits(mixobj2, "mix")) {
                if (sum(dim(mixobj2$mixdata) != dim(mixobj1$mixdata)) == 
                  0) {
                  if (sum(mixobj2$mixdata[, 2] != mixobj1$mixdata[, 
                    2]) == 0) {
                    nan1 <- sum(sapply(mixobj1$se, is.nan))
                    nan2 <- sum(sapply(mixobj2$se, is.nan))
                    if (nan1 != 0) {
                      warning(paste("There are ", nan1, "NaN values in se for model 1, ANOVA may be invalid"))
                      mixobj1$se <- apply(mixobj1$se, c(1, 2), 
                        function(x) if (is.nan(x)) 
                          -1
                        else x)
                    }
                    if (nan2 != 0) {
                      warning(paste("There are ", nan2, "NaN values in se for model 2, ANOVA may be invalid"))
                      mixobj2$se <- apply(mixobj2$se, c(1, 2), 
                        function(x) if (is.nan(x)) 
                          -1
                        else x)
                    }
                    if (sum(is.na(mixobj2$se) & !is.na(mixobj1$se)) == 
                      0 | sum(is.na(mixobj1$se) & !is.na(mixobj2$se)) == 
                      0) {
                      ddf <- max(sum(!is.na(mixobj2$se) & is.na(mixobj1$se)), 
                        sum(!is.na(mixobj1$se) & is.na(mixobj2$se)))
                      if (mixobj1$constraint$conpi == "PFX" | 
                        mixobj2$constraint$conpi == "PFX") 
                        ddf <- ddf - 1
                      if (ddf > 0) {
                        hypo <- TRUE
                        if (mixobj1$df > mixobj2$df) 
                          mixobj <- mixobj2
                        else mixobj <- mixobj1
                        if (ddf != max(mixobj1$df - mixobj2$df, 
                          mixobj2$df - mixobj1$df)) 
                          warning("Don't have same empty bins, or might not be nested")
                        dchisq <- max(mixobj1$chisq - mixobj2$chisq, 
                          mixobj2$chisq - mixobj1$chisq)
                      }
                      else warning("The models are not nested")
                    }
                  }
                }
            }
        }
        if (hypo) 
            mixanova <- data.frame(Df = c(ddf, mixobj$df), Chisq = c(dchisq, 
                mixobj$chisq), `Pr(>Chisq)` = c(pchisq(dchisq, 
                ddf, low = FALSE), mixobj$P), row.names = c("Hypothesis", 
                "Residuals"), check.names = FALSE)
        else mixanova <- data.frame(Df = mixobj1$df, Chisq = mixobj1$chisq, 
            `Pr(>Chisq)` = mixobj1$P, row.names = "Residuals", 
            check.names = FALSE)
    }
    else if (!missing(mixobj2)) {
        if (inherits(mixobj2, "mix")) 
            mixanova <- data.frame(Df = mixobj2$df, Chisq = mixobj2$chisq, 
                `Pr(>Chisq)` = mixobj2$P, row.names = "Residuals", 
                check.names = FALSE)
    }
    if (is.null(mixanova)) {
        warning("The class of argument(s) is not 'mix'")
        NULL
    }
    else {
        cat("\n")
        print(structure(mixanova, heading = "Analysis of Variance Table\n", 
            class = c("anova", class(mixanova))))
        cat("\n")
    }
}
