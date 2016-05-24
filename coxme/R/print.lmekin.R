print.lmekin <- function(x, ...) {
    if (x$method=="ML")
        cat("Linear mixed-effects kinship model fit by maximum likelihood\n")
    else cat("Linear mixed-effects kinship model fit by REML\n")
    cat("  Data:", deparse(x$call$data), "\n")
    if(!is.null(x$call$subset)) {
        cat("  Subset:", deparse(x$call$subset), "\n")
        }
    cat("  Log-likelihood =", format(x$loglik), "\n")
    omit <- x$na.action
    if(length(omit))
        cat("  n=", x$n, " (", naprint(omit), ")\n", sep = "")
    else cat("  n=", x$n, "\n\n")

    cat ("\nModel: ", deparse(x$call$formula), "\n")

    fcoef <- x$coefficients$fixed
    if (length(fcoef) > 0)  { # Not a ~1 model
        se <- sqrt(diag(x$var))
        tmp <- cbind(fcoef, se, round(fcoef/se,2),
               signif(1 - pchisq((fcoef/ se)^2, 1), 2))
        dimnames(tmp) <- list(names(fcoef), c("Value", "Std Error", "z", "p")) 

        cat("Fixed coefficients\n")
        print(tmp)
        }

    cat("\nRandom effects\n")

    random <- x$vcoef
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
#    temp1[,1] <- temp1[,1] * x$sigma
#    temp1[,2] <- temp1[,2] * x$sigma^2
    if (maxcol==3) temp1[,3] <- temp1[,3]*x$sigma

    indx <- cumsum(c(1, nrow))   # starting row of each effect
    temp3 <- rep("", nrow(temp1))
    temp3[indx[-length(indx)]] <- names(random)
    xname <- unlist(lapply(random, 
                 function(x) if (is.matrix(x)) dimnames(x)[[1]] else names(x)))
    temp <- cbind(temp3, xname, 
                  ifelse(is.na(temp1), "", format(temp1)))
    if (maxcol == 2)
        temp4 <- c("Group", "Variable", "Std Dev", "Variance")
    else 
        temp4 <- c("Group","Variable", "Std Dev", "Variance", "Corr", 
                   rep("", maxcol-3))
    dimnames(temp) <- list(rep("", nrow(temp)), temp4)
    print(temp, quote=F)
    cat("Residual error=", format(x$sigma), "\n")
    invisible(x)
    }
