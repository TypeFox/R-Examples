print.aftreg <- function(x, digits=max(options()$digits - 4, 3), ...){

    if (!is.null(cl<- x$call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
    }
    if (!is.null(x$fail)) {
        cat(" aftreg failed.\n")
        return()
    }
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    if (x$pfixed){

        n.slsh <- 1

    }else{
        n.slsh <- 2 * x$n.strata

    }
    coef <- x$coefficients
    se <- sqrt(diag(x$var))

    wald.p <- formatC(1 - pchisq((coef/ se)^2, 1),
                      digits = 3,
                      width = 9, format = "f")
    if (is.null(coef) || is.null(se))
        stop("Input is not valid")
#####################################
    if (x$param == "lifeAcc"){
        cat("Covariate          W.mean      Coef Time-Accn  se(Coef)    Wald p\n")
    }else{
        cat("Covariate          W.mean      Coef Life-Expn  se(Coef)    Wald p\n")
    }
    e.coef <- formatC(exp(coef), width = 9, digits = 3, format = "f")
    coef <- formatC(coef, width = 9, digits = 3, format = "f")
    se <- formatC(se, width = 9, digits = 3, format = "f")

    ett <-  formatC(1, width = 9, digits = 0, format = "f")
    noll <-  formatC(0, width = 5, digits = 0, format = "f")

    factors <- attr(x$terms, "factors")
    resp <- attr(x$terms, "response")
    row.strata <- attr(x$terms, "specials")$strata
    if (!is.null(row.strata))
        col.strata <- which(factors[row.strata, ] == 1)
    else col.strata <- NULL
    if (!is.null(x$covars)){
        if (!is.null(col.strata)){
            factors <-
              attr(x$terms, "factors")[-c(resp, row.strata), -col.strata,
                                       drop = FALSE]
        }else{
            factors <-
              attr(x$terms, "factors")[-c(resp, row.strata), ,
                                       drop = FALSE]
        }

        covar.names <- c(x$covars,
                         names(coef)[(length(coef)-n.slsh + 1):
                                       length(coef)])
        term.names <- colnames(factors)

        isF <- x$isF
    }
    ord <- attr(x$terms, "order")
    if (!is.null(col.strata)) ord <- ord[-col.strata]

    index <- 0

    if (!is.null(x$covars)){
        n.rows <- length(term.names)
        for (term.no in 1:n.rows){

            if (ord[term.no] == 1){
                covar.no <- which(factors[, term.no] == 1)

                if (isF[covar.no]){ ## Factors:
                    cat(covar.names[covar.no], "\n")

                    no.lev <- length(x$levels[[covar.no]])
                    cat(formatC(x$levels[[covar.no]][1], width = 16, flag =
                                "+"),
                        formatC(x$w.means[[covar.no]][1],
                                width = 8, digits = 3, format = "f"),
                        noll,
                        ett,
                        "          (reference)\n")
                    for (lev in 2:no.lev){

                        index <- index + 1
                        cat(formatC(x$levels[[covar.no]][lev], width = 16,
                                    flag = "+"),
                            formatC(x$w.means[[covar.no]][lev],
                                    width = 8, digits = 3, format = "f"),
                            coef[index],
                            e.coef[index],
                            se[index],
                            #formatC(" ", width = 1),
                            #formatC(" ", width = 9),
                            formatC(wald.p[index],
                                    digits = 3,
                                    width = digits + 2,
                                    format = "f"),
                            ##signif(1 - pchisq((coef/ se)^2, 1), digits - 1),
                            "\n")
                    }
                }else{ ## Covariates:
                    index <- index + 1
                    cat(formatC(covar.names[covar.no], width = 16, flag = "-"),
                        formatC(x$w.means[[covar.no]],
                                width = 8, digits = 3, format = "f"),
                        coef[index],
                        e.coef[index],
                                        #exp(coef[index]),
                        se[index],
                        #formatC(" ", width = 1),
                        formatC(wald.p[index],
                                digits = 3,
                                width = digits + 2,
                                format = "f"),
                        ##signif(1 - pchisq((coef/ se)^2, 1), digits - 1),
                        "\n")
                }
            }else if (ord[term.no] > 1){ ## Interactions:
                cat(format(term.names[term.no], width = 16), "\n")
                niv <- numeric(ord[term.no])
                covar.no <- which(factors[, term.no] == 1)

                for (i in 1:ord[term.no]){
                    if (isF[covar.no[i]]){
                        niv[i] <- length(x$levels[[covar.no[i]]]) - 1
                    }else{
                        niv[i] <- 1
                    }
                }
                stt <- index + 1
                for (index in stt:(stt + prod(niv) - 1)){
                    vn <- sub(covar.names[covar.no[1]], "", names(coef)[index])
                    for (i in 1:ord[term.no]){
                        vn <- sub(covar.names[covar.no[i]], "", vn)
                    }

                    cat(formatC(" ", width = 2),
                        formatC(substring(vn, 1, 22), width = 22, flag = "-"),
                        ## format(" ", 8),
                        coef[index],
                        e.coef[index],
                        se[index],
                        #formatC(" ", width = 1),
                        formatC(wald.p[index],
                                digits = 3,
                                width = digits + 2,
                                format = "f"),
                        ##signif(1 - pchisq((coef[index]/ se[index])^2, 1), digits - 1),
                        "\n")
                }
            }

        }
        cat("\n")
    }
    cat("Baseline parameters:\n")
    for (i in 1:n.slsh){
        jup <- length(coef)
        ss.names <- names(coef[(jup - n.slsh + 1):jup])
        index <- index + 1
        ## covar.no <- covar.no + 1
        cat(formatC(ss.names[i], width = 16, flag = "-"),
            formatC(" ",
                    width = 8, digits = 3, format = "c"),
            coef[index],
            formatC(" ", width = 9), # Changed 2.4-3
            
            ##e.coef[index],
                                        #exp(coef[index]),
            se[index],
            #formatC(" ", width = 1),
            formatC(wald.p[index],
                    digits = 3,
                    width = digits + 2,
                    format = "f"),
            ##signif(1 - pchisq((coef/ se)^2, 1), digits - 1),
            "\n")
    }
    cat("Baseline life expectancy: ", x$baselineMean, "\n")
#####################################
    if(FALSE){
        tmp <- cbind(coef, exp(coef), se,
                     signif(1 - pchisq((coef/ se)^2, 1), digits - 1))
        dimnames(tmp) <- list(names(coef), c("coef", "rel. risk",
                                             "se(coef)", "p"))

        cat("\n")
        prmatrix(tmp)
    }
    logtest <- -2 * (x$loglik[1] - x$loglik[2])
    if (is.null(x$df)) df <- sum(!is.na(coef)) - n.slsh
    else  df <- round(sum(x$df), 2)
    cat("\n")
    if (x$pfixed) {
        cat(" Shape is fixed at ", x$shape, "\n\n")
    }
    cat(formatC("Events", width = 25, flag = "-"), x$events, "\n")
    cat(formatC("Total time at risk", width = 25, flag = "-"),
        formatC(x$ttr, digits = 5, format = "fg"), "\n")
    cat(formatC("Max. log. likelihood", width = 25, flag = "-"),
        formatC(x$loglik[2], digits = 5, format = "fg"), "\n")
    if (df > 0.5){
        cat(formatC("LR test statistic", width = 25, flag = "-"),
            format(round(logtest, 2)), "\n")
        cat(formatC("Degrees of freedom", width = 25, flag = "-"),
            formatC(df, digits = 0, format = "f"), "\n")

        cat(formatC("Overall p-value", width = 25, flag = "-"),
          format.pval(1 - pchisq(logtest, df), digits = 6, "\n"))
    }
    cat("\n")
    if (length(x$icc))
      cat("   number of clusters=", x$icc[1],
          "    ICC=", format(x$icc[2:3]), "\n")
    invisible(x)
 }
