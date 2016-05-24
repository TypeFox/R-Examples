
`summary.p3state`<-function (object, model = NULL, covmat = NULL, estimate = NULL, 
    time1 = NULL, time2 = NULL, ...) 
{
    if (missing(object)) 
        stop("Argument 'object' is missing with no default")
    if (!inherits(object, "p3state")) 
        stop("'object' must be of class 'p3state'")
    if (missing(covmat)) 
        covmat <- FALSE
    if (missing(model)) 
        model <- "NONE"
    if (missing(time1)) 
        time1 <- 0
    if (missing(estimate) & missing(time2)) 
        estimate <- FALSE
    if (missing(time2)) 
        time2 <- max(object$datafr[, 1])
    if (missing(estimate)) 
        estimate <- TRUE
    if (time1 > time2) 
        stop("Argument 'time1' cannot be greater then 'time2'")
    if (time1 < 0 | time2 < 0) 
        stop("'time1' and 'time2' must be positive")
    if (object$descriptives[2] == 0) {
        cat("Progressive three-state model", "\n")
    }
    if (object$descriptives[2] > 0) {
        cat("Illness-death model", "\n")
    }
    cat("", "\n")
    cat("Number of individuals experiencing the intermediate event: ", 
        object$descriptives[1], "\n")
    cat("Number of events for the direct transition from state 1 to state 3: ", 
        object$descriptives[2], "\n")
    cat("Number of individuals remaining in state 1: ", object$descriptives[3], 
        "\n")
    cat("Number of events on transition from state 2: ", object$descriptives[4], 
        "\n")
    cat("Number of censored observations on transition from state 2: ", 
        object$descriptives[5], "\n")
    cat("", "\n")
    if (estimate == TRUE) {
	transprob<-pLIDA(object$datafr, time1, time2,tp="all")
       
        cat("The estimate of the transition probability P11(", 
            time1, ",", time2, ") is ", transprob[[1]], "\n")
        cat("The estimate of the transition probability P12(", 
            time1, ",", time2, ") is ", transprob[[2]], "\n")
	cat("The estimate of the transition probability P13(", 					#Changed
            time1, ",", time2, ") is ", 1 - transprob[[1]] - transprob[[2]], "\n")		#Changed
        cat("The estimate of the transition probability P22(", 
            time1, ",", time2, ") is ", transprob[[3]], "\n")
        cat("The estimate of the transition probability P23(", 
            time1, ",", time2, ") is ", 1 - transprob[[3]], "\n")
        if (object$descriptives[2] == 0) {
            biv <- Biv(object$datafr, time1, time2)
            cat("The estimate of the bivariate distribution function F12(", 
                time1, ",", time2, ") is ", biv, "\n")
            p7 <- which(object$datafr[, 3] <= time2 & object$datafr[, 
                2] * object$datafr[, 5] == 1)
            marg <- sum(object$datafr[p7, 6])
            cat("The estimate of the marginal distribution function of the second gap time, F2(", 
                time2, ") is ", marg, "\n")
        }
    }
    if (model == "TDCM") {
        cat("", "\n")
        cat("  ***** TIME-DEPENDENT COX REGRESSION MODEL *****  ", 
            "\n")
        cat("n= ", summary(object$tdcm)$n, "\n")
        print(summary(object$tdcm)$coef)
        cat(" ", "\n")
        print(summary(object$tdcm)$conf.int)
        cat(" ", "\n")
        cat("Likelihood ratio test= ", summary(object$tdcm)$logtest[[1]], 
            "on ", summary(object$tdcm)$logtest[[2]], " df, p=", 
            summary(object$tdcm)$logtest[[3]], "\n")
        cat(" ", "\n")
        cat("-2*Log-likelihood=", -2 * object$tdcm$loglik[2], 
            "\n")
        if (covmat == TRUE) {
            cat("*** variance-covariance matrix ***", "\n")
            print(object$tdcm$var)
            cat("\n")
        }
    }
    if (model == "CMM" | model == "CSMM") {
        cat("", "\n")
        if (model == "CMM") {
            cat("*********************** COX MARKOV MODEL ***********************", 
                "\n")
        }
        if (model == "CSMM") {
            cat("********************* COX SEMI-MARKOV MODEL *********************", 
                "\n")
        }
        cat("", "\n")
        if (object$descriptives[2] > 0) {
            cat("    *************** FROM STATE 1 TO STATE 3 ****************", 
                "\n")
            cat("", "\n")
            if (object$descriptives[2] < 5) 
                cat("Warning: there are few events on this transition", 
                  "\n")
            cat("n= ", summary(object$msm13)$n, "\n")
            print(summary(object$msm13)$coef)
            cat(" ", "\n")
            print(summary(object$msm13)$conf.int)
            cat(" ", "\n")
            cat("Likelihood ratio test= ", summary(object$msm13)$logtest[[1]], 
                "on ", summary(object$msm13)$logtest[[2]], " df, p=", 
                summary(object$msm13)$logtest[[3]], "\n")
            cat(" ", "\n")
            cat("-2*Log-likelihood=", -2 * object$msm13$loglik[2], 
                "\n")
            cat("\n")
            if (covmat == TRUE) {
                cat("*** variance-covariance matrix ***", "\n")
                print(object$msm13$var)
                cat("\n")
            }
        }
        cat("    *************** FROM STATE 1 TO STATE 2 ****************", 
            "\n")
        cat("", "\n")
        if (object$descriptives[1] < 5) 
            cat("Warning: there are few events on this transition", 
                "\n")
        cat("n= ", summary(object$msm12)$n, "\n")
        print(summary(object$msm12)$coef)
        cat(" ", "\n")
        print(summary(object$msm12)$conf.int)
        cat(" ", "\n")
        cat("Likelihood ratio test= ", summary(object$msm12)$logtest[[1]], 
            "on ", summary(object$msm12)$logtest[[2]], " df, p=", 
            summary(object$msm12)$logtest[[3]], "\n")
        cat(" ", "\n")
        cat("-2*Log-likelihood=", -2 * object$msm12$loglik[2], 
            "\n")
        cat("\n")
        if (covmat == TRUE) {
            cat("*** variance-covariance matrix ***", "\n")
            print(object$msm12$var)
            cat("\n")
        }
        if (model == "CMM") {
            cat("    *************** FROM STATE 2 TO STATE 3 ****************", 
                "\n")
            cat("", "\n")
            if (object$descriptives[4] < 5) 
                cat("Warning: there are few events on this transition", 
                  "\n")
            cat("n= ", summary(object$cmm23)$n, "\n")
            print(summary(object$cmm23)$coef)
            cat(" ", "\n")
            print(summary(object$cmm23)$conf.int)
            cat(" ", "\n")
            cat("Likelihood ratio test= ", summary(object$cmm23)$logtest[[1]], 
                "on ", summary(object$cmm23)$logtest[[2]], " df, p=", 
                summary(object$cmm23)$logtest[[3]], "\n")
            cat(" ", "\n")
            cat("-2*Log-likelihood=", -2 * object$cmm23$loglik[2], 
                "\n")
            cat("\n")
            if (covmat == TRUE) {
                cat("*** variance-covariance matrix ***", "\n")
                print(object$cmm23$var)
                cat("\n")
            }
        }
        if (model == "CSMM") {
            cat("    *************** FROM STATE 2 TO STATE 3 ****************", 
                "\n")
            cat("", "\n")
            if (object$descriptives[4] < 5) 
                cat("Warning: there are few events on this transition", 
                  "\n")
            cat("n= ", summary(object$csmm23)$n, "\n")
            print(summary(object$csmm23)$coef)
            cat(" ", "\n")
            print(summary(object$csmm23)$conf.int)
            cat(" ", "\n")
            cat("Likelihood ratio test= ", summary(object$csmm23)$logtest[[1]], 
                "on ", summary(object$csmm23)$logtest[[2]], " df, p=", 
                summary(object$csmm23)$logtest[[3]], "\n")
            cat(" ", "\n")
            cat("-2*Log-likelihood=", -2 * object$csmm23$loglik[2], 
                "\n")
            cat("\n")
            if (covmat == TRUE) {
                cat("*** variance-covariance matrix ***", "\n")
                print(object$csmm23$var)
                cat("\n")
            }
        }
        cat("Checking the Markov assumption:", "\n")
        cat("Testing if the time spent in state 1 (start) is important on transition from state 2 to state 3", 
            "\n")
        cat("\n")
        print(summary(object$tma)$coef)
        cat("\n")
        if (model == "CSMM") {
            if (summary(object$tma)$coef[1, 5] > 0.05) 
                cat("Warning: the p-value is ", summary(object$tma)$coef[1, 
                  5], "\n")
            if (summary(object$tma)$coef[1, 5] < 0.05) 
                cat("The p-value is ", summary(object$tma)$coef[1, 
                  5], "less than 5%", "\n")
            cat("\n")
        }
        if (model == "CMM") {
            if (summary(object$tma)$coef[1, 5] < 0.05) 
                cat("Warning: the p-value is ", summary(object$tma)$coef[1, 
                  5], "less than 5%", "\n")
            if (summary(object$tma)$coef[1, 5] > 0.05) 
                cat("The p-value is ", summary(object$tma)$coef[1, 
                  5], "\n")
            cat("\n")
        }
    }
}

