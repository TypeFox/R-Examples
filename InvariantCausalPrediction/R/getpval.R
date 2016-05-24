getpval <-
function(Y,X,IN,test="correlation",maxNoObs=200){

    
    linm <- glm( Y ~ X -1, family= if(is.factor(Y)) "binomial" else "gaussian", control=glm.control(maxit=10,epsilon=10^(-6)))
    coefficients <- coef(linm)[-1]
    coefficientsvar <- summary(linm)$coefficients[,2][-1]
    pred <- predict(linm, type="response")

    if(is.factor(Y)){
        resid <- ((as.numeric(Y)-1)-pred)
        pvalvec <- length(IN)
        for (ki in 1:length(IN)){
            pvalvec[ki] <- getpvalClassif(resid[IN[[ki]]],resid[-IN[[ki]]],test=test)
        }
        pval <- min(pvalvec)*(length(IN)-1)

    }else{
        useexact <- FALSE
        if( !is.function(test)){
            if( test=="exact") useexact <- TRUE
        }
        if(!useexact){
            
            K <- length(IN)
            
            resid <- residuals(linm)
            pvalvec <- numeric(length(IN))
            for (ki in 1:length(IN)){
                pvalvec[ki] <- pvalfunc( resid[IN[[ki]]], resid[-IN[[ki]]],test=test)
                if(!is.function(test)){
                    if(test=="correlation"){
                        nonzerosd <- which(apply(X,2,sd)>0)
                        if(length(nonzerosd)>0){
                            corpval <- numeric(length(nonzerosd))
                            for (k in 1:length(corpval)) corpval[k] <- cor.test(X[IN[[ki]],nonzerosd[k]], resid[IN[[ki]]])$p.value
                            pvalvec[ki] <- 2*min(pvalvec[ki], length(corpval)*min(corpval))
                        }
                    }
                }
            }
            pval <- min(pvalvec)*(length(IN)-1)
            
        }else{
            
            K <- length(IN)
            
            if(!is.matrix(X)) X <- as.matrix(X)
            n <- nrow(X)
            
            pvalvec <- numeric(length(IN))
            for (ki in 1:length(IN)){
                nk <- length(IN[[ki]])
                nko <- n-nk
                p <- ncol(X)
                linm <- lm.fit( X[-IN[[ki]], ,drop=FALSE] , Y[-IN[[ki]]]) ## fit a model on all other data
                pred <- as.numeric(X[IN[[ki]], ,drop=FALSE] %*% coefficients(linm))
                diff <- Y[IN[[ki]]] - pred
                
                selobs <-  if( nk>maxNoObs)  sample( IN[[ki]], maxNoObs) else IN[[ki]]
                if( nk>maxNoObs) diff <- diff[ IN[[ki]] %in% selobs]
                nk <- length(selobs)
                COV <- diag(length(diff)) + X[ selobs,] %*% solve(t(X[-IN[[ki]],])%*%X[-IN[[ki]],], t(X[ selobs,]))
                
                stat <- (t(diff)%*% solve(COV, diff)) / (nk * var(residuals(linm))*nko/(nko-p))
                pval <- 1-pf(stat, nk, n-nk - ncol(X))
                
                pvalvec[ki] <- pval
            }
        }
        
        pval <- min(pvalvec) * (length(IN))
    }
    pval <- min(1,pval)
    
    return(list(pval=pval,coefficients=coefficients,coefficientsvar=coefficientsvar))
}
