# datamat is a numeric data matrix (exclude the groups from the first column)
LinearModelFit <- function(datamat, 
    factormat=NULL,
    covariatemat=NULL,
    contrastmat=NULL,
    ruv2=TRUE,
    k=NULL, nc=NULL, # These are as defined in the Normalise function
    moderated=FALSE,
    padjmethod="BH",
    saveoutput=FALSE,
    outputname="results", ...)
{
    datamat <- as.matrix(datamat)
    if (mode(datamat)!="numeric" ) {
        stop("datamat must be a numerical data matrix")
    }
    Y <- t(datamat)
    
    if (ruv2){
        # If no k, get them to enter it
        if (is.null(k)) {
            stop("Please enter the number of unwanted variation factors")
        }
        # If there is no nc, get them to enter it
        if (is.null(nc)) {
            stop(
                paste("Please enter a vector with columns",
                    "corresponding to non-changing metabolites"
                )
            )
        }
        if (!is.null(covariatemat)){
            if (class(covariatemat)!="matrix") {
                stop("covariatemat must be a matrix")
            }
            z <- covariatemat
            rzY <- t(Y) - z %*% solve(t(z) %*% z) %*% t(z) %*% t(Y) 
            factormat <- factormat - z %*% solve(t(z) %*% z) %*% t(z) %*%
                factormat
        } else {
            rzY <- t(Y)
        }
        
        Wruv2 <- svd(rzY[, nc] %*% t(rzY[, nc]))$u[, 1:k, drop = FALSE]
        
        colnames(Wruv2) <- paste("k", c(1:k), sep="")
        designmat <- cbind(factormat, covariatemat, Wruv2)
        if (!is.null(contrastmat)) {
            contrastmat <- rbind(
                contrastmat, matrix(0, nrow=k, ncol=ncol(contrastmat))
            )
            rownames(contrastmat) <- colnames(designmat)
        }
    } else {
        designmat <- cbind(factormat, covariatemat)
    }
    #Linear model fit
    fit <- lmFit(Y, design=designmat, ...)
    #Contrast fits if required
    if (!is.null(contrastmat)) {
        fit <- contrasts.fit(fit=fit, contrasts=contrastmat, ...)
    }
    #eBayes fit
    ebfit <- eBayes(fit, ...)
    ebfit$metabolites <- ebfit$genes
    ebfit$genes <- ebfit$rank <- ebfit$assign <- NULL
    ebfit$qr <- ebfit$qraux <- ebfit$pivot <- ebfit$tol <- NULL
    ebfit$cov.coefficients <- ebfit$pivot <- ebfit$lods <- NULL


    #statistics
    beta <- ebfit$coeff
    if (moderated) {
        Fstat <- ebfit$F
        Fpval <- ebfit$F.p.value
        tstat <- ebfit$t
        tpval <- ebfit$p.value
        se <- beta / tstat
    } else {
        tstat <- sweep(
            data.matrix(ebfit$coef / ebfit$stdev.unscaled), 1, ebfit$sigma, "/"
        )
        tpval <- 2 * pt(-abs(tstat), df=ebfit$df.residual)
        ebfit$t <- tstat
        df <- ebfit$df.residual
        fstat <- classifyTestsF(ebfit, df=df, fstat.only=TRUE)
        Fstat <- as.vector(fstat)
        df1 <- attr(fstat, "df1")
        df2 <- attr(fstat, "df2")
        if (df2[1] > 1e+06) {
            Fpval <- pchisq(df1 * Fstat, df1, lower.tail=FALSE)
        } else {
            Fpval<- pf(Fstat, df1, df2, lower.tail=FALSE)
        }
        se <- beta / tstat
        
        ebfit$F <- Fstat
        ebfit$F.p.value <- Fpval
        ebfit$p.value <- tpval
        ebfit$df.total <- ebfit$s2.post <- ebfit$stdev.unscaled <- NULL
        ebfit$var.prior <-ebfit$proportion <- ebfit$s2.prior <- NULL
        ebfit$df.prior <- NULL
        ebfit$std.error <- se
        ebfit$t <- tstat
        ebfit$df <- df
    }
    
    # ruvmat
    if (ruv2) {
        ebfit$uvmat <- Wruv2
    }
    # adjusted p for t p value
    tpadj <- matrix(NA,ncol=ncol(tstat),nrow=nrow(tstat))
    colnames(tpadj) <- colnames(tpval)
    row.names(tpadj) <- row.names(tpval)
    for (j in 1:ncol(tpadj)) {
        tpadj[,j] <- p.adjust(tpval[,j], method=padjmethod, 
            n=length(na.omit(tpval[,j]))
        )
    }
    ebfit$adj.p.value <- tpadj 
    
    mat <- data.frame(Fstat,
        Fpval,
        p.adjust(Fpval, method=padjmethod, n=length(na.omit(Fpval))),
        beta, tstat, se, tpval, tpadj
    )
    
    colnames(mat) <- c("F stat",
        "F p value",
        "Adjusted F p value",
        paste("coeff", colnames(ebfit$coeff)),
        paste("t stat", colnames(ebfit$t)),
        paste("std error", colnames(ebfit$coeff)),
        paste("t p value", colnames(ebfit$t)),
        paste("Adjusted t p value", colnames(ebfit$t))
    )
    
    if (saveoutput) {
        write.csv(mat,paste(c(outputname,".csv"),collapse=""))
    }
    
    return(fit=ebfit)
}
