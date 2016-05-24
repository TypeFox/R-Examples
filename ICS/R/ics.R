### function to constract ICS or ICA based on two scatter matrices
### returns an object of S4 class ics
###

`ics` <-
    function (X, S1 = cov, S2 = cov4, S1args = list(), S2args = list(), 
    stdB = "Z", stdKurt=TRUE, na.action = na.fail) 
    {
        X <- na.action(X)
        data.matrix <- as.matrix(X)
        
        if (stdB != "B" & stdB != "Z") 
            stop("'stdB' must be 'B' or 'Z'")
        p <- dim(X)[2]
        if (p < 2) 
            stop("'X' must be at least bivariate")
        if (!is.logical(stdKurt)) stop("'stdKurt must be logical'")
    
        S1name <- deparse(substitute(S1))
        S2name <- deparse(substitute(S2))
    
        if (is.matrix(S1))
        {
        if (!is.matrix(S2)) stop("'S1' and 'S2' must be both functions or both square matrices")
        p1<-dim(S1)
        if (!isSymmetric(S1)) stop("'S1' must be a symmetric matrix")
        if (p1[1]!=p) stop("'S1' must be a square matrix that fits to the data")
        p2<-dim(S2)
        if (!isSymmetric(S2)) stop("'S2' must be a symmetric matrix")
        if (p1[1]!=p2[1]) stop("'S1' and 'S2' must have the same dimension")
    
        B1 <- solve(mat.sqrt(S1))
        X1 <- data.matrix %*% B1
        B2 <- B1 %*% S2 %*% B1
        }
    
        if (is.function(S1))
        {
        if (!is.function(S2)) stop("'S1' and 'S2' must be both functions or both symmetric matrices")
        S1 <- do.call("S1", c(list(X), S1args))
        B1 <- solve(mat.sqrt(S1))
        X1 <- data.matrix %*% B1
        S2 <- do.call("S2", c(list(X), S2args))
        B2 <- B1 %*% S2 %*% B1
        }
    
    
        B2.eigen <- eigen(B2)
        U2 <- B2.eigen$vectors
        DiagB2 <- B2.eigen$values
        if (stdKurt == TRUE) DiagB2 <- DiagB2/prod(DiagB2)^(1/p)
        X2 <- X1 %*% U2
        B <- crossprod(U2, B1)
        
        # choosing the signs of B
        
        if (stdB == "B") {
            row.signs <- apply(B, 1, .sign.max)
            row.norms <- sqrt(rowSums((B)^2))
            B.res <- sweep(B, 1, row.norms * row.signs, "/")
            Z <- as.data.frame(tcrossprod(data.matrix, B.res))
        }
        if (stdB == "Z") {
            Z1 <- tcrossprod(data.matrix, B)
            skewness <- colMeans(Z1) - apply(Z1, 2, median)
            skew.signs <- ifelse(skewness > 0, 1, -1)
            B.res <- sweep(B, 1, skew.signs, "*")
            Z <- as.data.frame(tcrossprod(data.matrix, B.res))
        }
        names(Z) <- paste(rep("IC", p), 1:p, sep = ".")
        if (is.null(colnames(X)) == TRUE) 
            names.X <- paste(rep("X", p), 1:p, sep = ".")
        else names.X <- colnames(X)
        res <- new("ics", gKurt = DiagB2, UnMix = B.res, S1 = S1, S2 = S2, S1name = S1name, 
            S2name = S2name, Scores = Z, DataNames = names.X, 
            StandardizeB = stdB, StandardizegKurt = stdKurt)
        return(res)
    }
