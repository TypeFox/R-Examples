########### R function: almeAux #################

# For obtaining df and conditional (on variance components) variance 
# of fixed and random parameter estimates and predictions in linear
# mixed models, modified compared to lmeAux to take into account 
# correlated residuals

# Last changed: 16 JUN 2006


"almeAux" <-
function (X, Z, G, RR, resid.var, block.inds = NA, indiv = 1, rho = 0) 
{
    n <- nrow(as.matrix(X))
    if (!is.null(Z)) {
        q <- ncol(as.matrix(Z))
        p <- ncol(as.matrix(X))
        if (is.list(block.inds) == 0) 
            block.inds <- list(1:(p + q))
        if (qr(G)$rank == ncol(G)) 
            G.tilde <- rbind(matrix(0, p, p + q), cbind(matrix(0, 
                q, p), resid.var * solve(G)))
        if (qr(G)$rank < ncol(G)) 
            G.tilde <- rbind(matrix(0, p, p + q), cbind(matrix(0, 
                q, p), diag(rep(1e+06, q))))
        C.mat <- cbind(X, Z)
    }
    if (is.null(Z)) 
        C.mat <- X
    CTC <- t(C.mat)%*%C.mat

#take into account given correlation matrix of the residuals.

    RR.inv <- NULL
    if(!is.null(RR))
      {
        RR <- chol(RR)
        rr.inv <- backsolve(RR,diag(rep(1,ncol(RR))))
        RR.inv <- rr.inv%*%t(rr.inv)
        CTC <- t(C.mat) %*% RR.inv%*%C.mat
       }
    CTRC <- acompute.CTRinvC(X, Z, RR.inv,indiv, rho)
    if (is.null(Z)) 
        G.tilde <- matrix(0, nrow(CTRC), ncol(CTRC))
    Ridge <- CTRC + G.tilde
    Ridge <- (Ridge + t(Ridge))/2
    R.ridge <- chol(Ridge)
    R.ridge.inv <- backsolve(R.ridge, diag(rep(1, nrow(R.ridge))))
    ridge.inv <- R.ridge.inv %*% t(R.ridge.inv)
    df.mat <- ridge.inv %*% CTC
    cov.mat <- resid.var * ridge.inv
    df <- numeric()
    for (j in 1:length(block.inds)) {
        curr.inds <- block.inds[[j]]
        df[j] <- sum(diag(as.matrix(df.mat[curr.inds, curr.inds])))
    }
    df.fit <- sum(df)
    df.res <- n - 2 * df.fit + sum(diag(df.mat %*% df.mat))
    lmeAux.object <- list(cov.mat = cov.mat, df = df, block.inds = block.inds, 
        resid.var = resid.var, random.var = G, df.fit = df.fit, 
        df.res = df.res)
    return(lmeAux.object)
}
