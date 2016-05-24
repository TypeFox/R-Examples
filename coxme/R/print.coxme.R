print.coxme <- function(x, rcoef=FALSE, digits=options()$digits, ...) {
    cat("Cox mixed-effects model fit by maximum likelihood\n")
    if (!is.null(x$call$data)) 
        cat("  Data:", deparse(x$call$data))
    if(!is.null(x$call$subset)) {
        cat(";  Subset:", deparse(x$call$subset), "\n")
	}
    else cat("\n")

    beta <- x$coefficients
    nvar <- length(beta)
    nfrail<- nrow(x$var) - nvar

    omit <- x$na.action
    cat("  events, n = ", x$n[1], ', ', x$n[2], sep='')
    if(length(omit))
        cat(" (", naprint(omit), ")", sep = "")
    loglik <- x$loglik + c(0,0, x$penalty)
    temp <- matrix(loglik, nrow=1)
    cat("\n  Iterations=", x$iter, "\n")
    dimnames(temp) <- list("Log-likelihood", 
                           c("NULL", "Integrated", "Fitted"))
    print(temp)
    cat("\n")
    chi1 <- 2*diff(x$loglik[c(1,2)]) 

    
    chi1 <- 2*diff(loglik[1:2]) 
    chi2 <- 2*diff(loglik[c(1,3)])
    temp <- rbind(c(round(chi1,2), round(x$df[1],2),
                    signif(1- pchisq(chi1,x$df[1]),5),
                    round(chi1- 2*x$df[1],2),
                    round(chi1- log(x$n[1])*x$df[1],2)),
                  c(round(chi2,2), round(x$df[2],2),
                    signif(1- pchisq(chi2,x$df[2]),5),
                    round(chi2- 2*x$df[2],2),
                    round(chi2- log(x$n[1])*x$df[2],2)))
    dimnames(temp) <- list(c("Integrated loglik", " Penalized loglik"),
                           c("Chisq", "df", "p", "AIC", "BIC"))
    print(temp, quote=F, digits=digits)

    cat ("\nModel: ", deparse(x$call$formula), "\n")

    if (nvar > 0)  { # Not a ~1 model
        se <- sqrt(diag(x$var)[nfrail+1:nvar])
        tmp <- cbind(beta, exp(beta), se, round(beta/se,2),
               signif(1 - pchisq((beta/ se)^2, 1), 2))
        dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)",
            "se(coef)", "z", "p"))
        }
    if (rcoef) { # print the random coefs
        #next line unlists, trying to give good names to the coefs
        coef <- unlist(lapply(ranef(x), function(y) {
            if (is.matrix(y)) {
                z <- c(y)
                dd <- dimnames(y)
                names(z) <- c(outer(dd[[1]], dd[[2]], paste,sep=':'))
                z}
            else y
            }))
                                  
        se <- sqrt(diag(x$var)[1:nfrail])
        rtmp <- cbind(coef, exp(coef), se)
        dimnames(rtmp) <- list(names(coef), c("coef", "exp(coef)",
            "Penalized se"))
        }

    if (nvar>0 && rcoef) {
        cat("Fixed and penalized coefficients\n")
        print(rbind(tmp, cbind(rtmp,NA,NA)), na.print='', digits=digits)
        }
    else if (rcoef) {
        cat("Penalized coefficients\n")
        print(rtmp, digits=digits)
        }
    else if (nvar>0) {
        cat("Fixed coefficients\n")
        print(tmp, digits=digits)
        }

    cat("\nRandom effects\n")

    random <- VarCorr(x)
    nrow <-  sapply(random, 
                    function(x) if (is.matrix(x)) nrow(x) else length(x))
    maxcol <-max(sapply(random,
                        function(x) if (is.matrix(x)) 1+ncol(x) else 2))
    temp1 <- matrix(NA, nrow=sum(nrow), ncol=maxcol)
    indx <- 0
    for (term in  random) {
        if (is.matrix(term)) {
            k <- nrow(term)
            nc <- ncol(term)  #assume nc > nr (only cases I know so far)
            for (j in 1:k) {
                temp1[j+indx, 1] <- sqrt(term[j,j])
                temp1[j+indx, 2] <- term[j,j]
                if (nc>j) {
                    indx2 <- (j+1):nc
                    temp1[j+indx, 1+ indx2] <- term[j, indx2]
                    }
                }
            }
        else {
            k <- length(term)
            temp1[1:k + indx,1] <- sqrt(term)
            temp1[1:k + indx,2] <- term
            }
        indx <- indx + k
        }
        
    indx <- cumsum(c(1, nrow))   # starting row of each effect
    temp3 <- rep("", nrow(temp1))
    temp3[indx[-length(indx)]] <- names(random)
    xname <- unlist(lapply(random, 
                  function(x) if (is.matrix(x)) dimnames(x)[[1]] else names(x)))
    temp <- cbind(temp3, xname, ifelse(is.na(temp1), "", 
                                       format(temp1, digits=digits)))
    if (maxcol == 2)
        temp4 <- c("Group", "Variable", "Std Dev", "Variance")
    else 
        temp4 <- c("Group","Variable", "Std Dev", "Variance", "Corr", 
                   rep("", maxcol-3))
    dimnames(temp) <- list(rep("", nrow(temp)), temp4)

    print(temp, quote=F)
    invisible(x)
    }

summary.coxme <- function(object, ...)
    print.coxme(object, ...)
