anova.mvl1lm <- function(object, object2=NULL, test = "Score", ...)
        {
        if (is.null(object2)){
            res <- list(models = object$call, statistic = object$statistic , parameter = object$parameter, p.value = object$p.value, method = "Score type test that all coefficients are 0:")
            class(res) <- "anovamvl1lm"
            return(res)
            }
        
        T.scores <- object$scores
        if (T.scores != object2$scores) stop("'object' and 'object2' must be fitted using the same scores 'x'")
        
        n <- dim(object$residuals)[1]
        p <- dim(object$residuals)[2]
        if (n != dim(object2$residuals)[1]) stop("'object' and 'object2' must be fitted using the same number of observations")
        
        test <- match.arg(test, c("Score", "Wald"))
        
        models <- c(object$call, object2$call)
        
        switch(test,
                "Score"={
                        if (any(c(is.null(object$x),is.null(object2$x)))) stop("'object' and 'object2' must contain the design matrix 'x'") 
                        
                        M.full <- object$x
                        M.nested <- object2$x
        
                        M.full.names <- colnames(M.full)
                        M.nest.names <- colnames(M.nested)
                        if (!all(M.nest.names %in% M.full.names)) stop("'object2' is not nested in 'object'") 
                        
        
                        X.2 <- subset(M.full,select = M.full.names[!(M.full.names %in% M.nest.names)]) 
                        X.1 <- subset(M.full,select = M.full.names[(M.full.names %in% M.nest.names)]) 
                        Y.hat <- object2$residuals
                        ch.XX.1 <- chol(crossprod(X.1))
                        P.X.1 <- X.1 %*% backsolve(ch.XX.1, forwardsolve(ch.XX.1, t(X.1), upper.tri=TRUE, transpose=TRUE)) #X.1 %*% solve(crossprod(X.1)) %*% t(X.1)
                        X.2.hat <- crossprod((diag(1,n) - P.X.1), X.2)
                        ch.XX.2.hat <- chol(crossprod(X.2.hat))
                        P.X.2.hat <- X.2.hat %*% backsolve(ch.XX.2.hat, forwardsolve(ch.XX.2.hat, t(X.2.hat), upper.tri=TRUE, transpose=TRUE)) #X.2.hat %*% solve(crossprod(X.2.hat)) %*% t(X.2.hat)
                        q.2 <- dim(X.2)[2]
                        method <- "Score type test that coefficients not in the restricted model are 0:"
                        switch(T.scores,
                          "identity" = {
                              # Q.2 <- n * sum(diag(crossprod(Y.hat,P.X.2.hat) %*% Y.hat %*% solve(crossprod(Y.hat))))
                               Q.2 <- n * sum(diag(crossprod(Y.hat,P.X.2.hat) %*% tcrossprod(Y.hat,syminv(crossprod(Y.hat)))))
                              dfs <- p*q.2
                              p.value <- 1 - pchisq(Q.2, df = dfs)
                              },
                          "sign" = {
                              switch(object$stand,
                                    "outer"={
                                            U.hat <- spatial.sign(Y.hat, center=FALSE, shape=FALSE)
                                            # Q.2 <- n * sum(diag(crossprod(U.hat,P.X.2.hat) %*% U.hat %*% solve(crossprod(U.hat))))
                                            Q.2 <- n * sum(diag(crossprod(U.hat,P.X.2.hat) %*% tcrossprod(U.hat, syminv(crossprod(U.hat)))))
                                            dfs <- p*q.2
                                            p.value <- 1 - pchisq(Q.2, df = dfs)
                                            },
                                    "inner"={
                                            U.hat <- spatial.sign(Y.hat, center=FALSE, shape=object2$S.mat)
                                            # Q.2 <- n * sum(diag(crossprod(U.hat,P.X.2.hat) %*% U.hat %*% solve(crossprod(U.hat))))
                                            Q.2 <- n * sum(diag(crossprod(U.hat,P.X.2.hat) %*% tcrossprod(U.hat, syminv(crossprod(U.hat)))))
                                            dfs <- p*q.2
                                            p.value <- 1 - pchisq(Q.2, df = dfs)
                                    }
                              )
                              },
                        
                          "rank" = {       
                              if (object$IntC != object2$IntC) stop("for rank regression anova cannot test the intercept term")                 
                              switch(object$stand,
                                    "outer"={
                                            R.hat <- spatial.rank(Y.hat, center=FALSE, shape=FALSE)
                                            # Q.2 <- n * sum(diag(crossprod(R.hat,P.X.2.hat) %*% R.hat %*% solve(crossprod(R.hat))))
                                            Q.2 <- n * sum(diag(crossprod(R.hat,P.X.2.hat) %*% tcrossprod(R.hat, syminv(crossprod(R.hat)))))
                                            dfs <- p*q.2
                                            p.value <- 1 - pchisq(Q.2, df = dfs)
                                            },
                                    "inner"={
                                            R.hat <- spatial.rank(Y.hat, center=FALSE, shape=object2$S.mat)
                                            # Q.2 <- n * sum(diag(crossprod(R.hat,P.X.2.hat) %*% R.hat %*% solve(crossprod(R.hat))))
                                            Q.2 <- n * sum(diag(crossprod(R.hat,P.X.2.hat) %*% tcrossprod(R.hat, syminv(crossprod(R.hat)))))
                                            dfs <- p*q.2
                                            p.value <- 1 - pchisq(Q.2, df = dfs)
                                    
                                    }
                              )
                              }
                              )
                              }
                ,
                "Wald"={
                        vcov.full <- object$vcov
                        vcov.nest <- object2$vcov
        
                        vcov.full.names <- colnames(vcov.full)
                        vcov.nest.names <- colnames(vcov.nest)
                        coef.full <- object$coef
                        coef.nest <- object2$coef
                        coef.full.names <- rownames(coef.full)
                        coef.nest.names <- rownames(coef.nest)
                        ind.C <- which( !(coef.full.names %in% coef.nest.names))
                        coef.diff <- coef.full[ind.C, ]
                        beta.diff <- as.vector(coef.diff)
                        method <- "Wald type test that coefficients not in the restricted model are 0:"
        
                        #p.Y <- dim(coef.full)[2]
                        #q.betaF <- dim(coef.full)[1]
                        #q.betaA <- dim(coef.diff)[1]
                        #ind.V <- rep(ind.C,p.Y) + (rep(1:q.betaA, each=p.Y)-1)*q.betaF
        
                        #if (!all(vcov.nest.names %in% vcov.full.names)) stop("'object2' is not nested in 'object'") 
                        ind.V <- which( !(vcov.full.names %in% vcov.nest.names))
                        vcov.diff <- vcov.full[ind.V, ind.V]
                        # Q.2 <- as.numeric(beta.diff %*% solve(vcov.diff) %*% beta.diff)
                        ch.vcov.diff <- chol(vcov.diff)
                        Q.2 <- as.numeric(beta.diff %*% backsolve(ch.vcov.diff, forwardsolve(ch.vcov.diff, beta.diff, upper.tri=TRUE, transpose=TRUE)))
                        dfs <- length(beta.diff)
                        p.value <- 1 - pchisq(Q.2, df = dfs)
                        }
                        
               )
        
        res <- list(models = models, method = method, statistic = Q.2 , parameter = dfs, p.value = p.value)
        class(res) <- "anovamvl1lm"
        return(res)
        }
    
