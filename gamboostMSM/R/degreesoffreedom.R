degreesoffreedom <- function(m, statusinfo = TRUE){
    P <- length(m$baselearner)
    ms <- m$mstop()
    mx <- m$xselect()
    N <- nrow(as.matrix(extract(m$baselearner[[1]], what="design", asmatrix=T, expand=T)))
    H <- vector("list", P)
    for(j in sort(unique(mx))){
        X <- extract(m$baselearner[[j]], what="design", asmatrix=T, expand=T)
        K <- extract(m$baselearner[[j]], what="penalty", asmatrix=T)
        lambda <- as.numeric(extract(m$baselearner[[j]]$dpp(rep(1, N)), what = "lambda"))
        H[[j]] <- diag(N) - nu*(X%*%solve((t(X)%*%X)+lambda*K)%*%t(X))
    }
    DF <- rep(0, ms)
    nu <- m$control$nu
    backpart <- diag(N)
    if(statusinfo){
        cat("  ")
        for(i in 1:ms){
            backpart <- backpart%*%H[[mx[i]]]
            DF[i] <- sum(diag(diag(N)-backpart))
            if(i %% 80 != 0){
                cat(".")
            }else{
                cat(paste(". df: ", round(DF[i], 2), sep=""))
                cat("\n  ")
            }
        }
    }else{
        for(i in 1:ms){
            backpart <- backpart%*%H[[mx[i]]]
            DF[i] <- sum(diag(diag(N)-backpart))
        }
    }
    return(DF)}
