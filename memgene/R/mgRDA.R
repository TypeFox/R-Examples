mgRDA <-
function(genD, vectorsMEM, perm=NULL, full=TRUE) {
    
    
    ## Internal function to do RDA on a genetic distance matrix
    .genRDA <- function(Y, X, full=FALSE) {
     
        n <- nrow(Y)
        p <- ncol(X)
        
        ## Double-centre genetic distance matrix
        row.wt = rep(1, nrow(Y))
        col.wt = rep(1, ncol(Y))
        st <- sum(col.wt)
        sr <- sum(row.wt)
        row.wt <- row.wt/sr
        col.wt <- col.wt/st
        Y <- -0.5*(Y*Y)
        row.mean <- apply(row.wt * Y, 2, sum)
        col.mean <- apply(col.wt *t(Y), 2, sum)
        col.mean <- col.mean - sum(row.mean * col.wt)
        Y <- sweep(Y, 2, row.mean)
        G <- t(sweep(t(Y), 2, col.mean))
        
        ## Regression
        H <- X %*% solve(t(X)%*%X) %*% t(X)
        I <- diag(n)
        predicted <- H%*%G%*%H
        res <- (I-H)%*%G%*%(I-H)
        MS_regression <- sum(diag(predicted))/p
        MS_residual <- sum(diag(res))/(n-p-1)
        F <- MS_regression/MS_residual
        MS_Total <- sum(diag(G))/n
        RsqAdj <- 1 - MS_residual/MS_Total
        
        if (full) {
            ## PCA on predicted values
            prC <- princomp(predicted)
            rdaPrC <- prC$scores
            rdaSdev <- prC$sdev
            dimnames(rdaPrC)[[2]] <- paste("MEMGENE", 1:ncol(rdaPrC), sep="")
            names(rdaSdev) <- paste("MEMGENE", 1:length(rdaSdev), sep="")
            
            return(list(RsqAdj=RsqAdj,
                    F=F,
                    resid=res,
                    pred=predicted,
                    memgene=rdaPrC,
                    sdev=rdaSdev))
        }
        else {
            return(list(RsqAdj=RsqAdj,
                    F=F,
                    resid=res,
                    pred=predicted))
        }
    }
    
    
    Y <- as.matrix(genD)
    X <- as.matrix(vectorsMEM)
    actual <- .genRDA(Y, X, full=full) 
    
    ## Perform permutation test if required
    if (!is.null(perm)) {
        n <- nrow(Y)
        RsqAdj <- actual$RsqAdj
        F <- actual$F
        Prob <- 1/perm;
        for (i in 1:(perm-1)) {
            permuted_rows <- sample(n, replace=FALSE)
            result <- .genRDA(Y, X[permuted_rows, ])
            if (result$F >= F) {
                Prob <- Prob + 1/perm
            }
        }
        actual$P <- Prob
        return(actual)
    
    }
    else {
        actual$P <- "Not tested"
        return(actual)
    }

}
