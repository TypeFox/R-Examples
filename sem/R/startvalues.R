# last modified 2013-01-24 by J. Fox

startvalues2 <- function(S, ram, debug=FALSE, tol=1E-6){
	ram <- ram[!(ram[, "parameter"] == 0 & ram[, "start value"] == 0), , drop=FALSE] # eliminate fixed parameters set to 0
    n <- nrow(S) 
    observed <- 1:n       
    m <- max(ram[,c(2,3)])            
    t <- max(ram[,4])   
    s <- sqrt(diag(S))
    R <- S/outer(s,s)
    latent<-(1:m)[-observed]
    par.posn <-  sapply(1:t, function(i) which(ram[,4] == i)[1])    
    one.head <- ram[,1] == 1
    start <- (ram[,5])[par.posn]
    A.pat <-matrix(FALSE, m, m)
    A.pat[ram[one.head, c(2,3), drop=FALSE]] <- TRUE
    P.pat <- C <- matrix(0, m, m)
    P.pat[ram[!one.head, c(2,3), drop=FALSE]] <- P.pat[ram[!one.head, c(3,2), drop=FALSE]] <- 1
    C[observed, observed] <- R
    for (l in latent) {
        indicators <- A.pat[observed, l]
        for (j in observed){
            C[j, l] <- C[l, j] <- if (!any(indicators)) runif(1, .3, .5)
                else {   
                        numerator <- sum(R[j, observed[indicators]])
                        denominator <- sqrt(sum(R[observed[indicators], observed[indicators]]))
                        numerator/denominator
                    }
            }
        }
    for (l in latent){
        for (k in latent){
            C[l, k] <- C[k,l] <- if (l==k) 1 else {
                                indicators.l <- A.pat[observed, l]
                                indicators.k <- A.pat[observed, k]
                                if ((!any(indicators.l)) | (!any(indicators.k))) runif(1, .3, .5) else {
                                    numerator <- sum(R[observed[indicators.l], observed[indicators.k]])
                                    denominator <- sqrt( sum(R[observed[indicators.l], observed[indicators.l]])
                                        * sum(R[observed[indicators.k], observed[indicators.k]]))
                                    numerator/denominator}
                                    }
            }
        }
    A <- matrix(0, m, m)
    for (j in 1:m){
        ind <- A.pat[j,]
        if (!any(ind)) next
        A[j, ind] <- solve(C[ind, ind]) %*% C[ind, j]
        }
    A[observed,] <- A[observed,]*matrix(s, n, m)
    A[,observed] <- A[,observed]/matrix(s, m, n, byrow=TRUE)
    C[observed,] <- C[observed,]*matrix(s, n, m)
    C[,observed] <- C[,observed]*matrix(s, m, n, byrow=TRUE)
    P <- (diag(m) - A) %*% C %*% t(diag(m) - A)
    P <- P.pat * P
    for (par in 1:t){
        if (!is.na(start[par])) next
        posn <- par.posn[par]
        if (ram[posn, 1] == 1) start[par] <- A[ram[posn, 2], ram[posn, 3]]
            else start[par] <- P[ram[posn, 2], ram[posn, 3]]
        if (abs(start[par]) < tol) start[par] <- tol
        }
    if (debug){
        cat('\nStart values:\n')
        print(start)
        cat('\n')
        }
    start
    }

startvalues <- function(S, ram, debug=FALSE, tol=1E-6){
    ram <- ram[!(ram[, "parameter"] == 0 & ram[, "start value"] == 0), , drop=FALSE] # eliminate fixed parameters set to 0
    n <- nrow(S) 
    observed <- 1:n       
    m <- max(ram[,c(2,3)])            
    t <- max(ram[,4])   
    s <- sqrt(diag(S))
    R <- S/outer(s,s)
    latent<-(1:m)[-observed]
    par.posn <-  sapply(1:t, function(i) which(ram[,4] == i)[1])    
    one.head <- ram[,1] == 1
    start <- (ram[,5])[par.posn]
    A.pat <-matrix(FALSE, m, m)
    A.pat[ram[one.head, c(2,3), drop=FALSE]] <- TRUE
    P.pat <- C <- matrix(0, m, m)
    P.pat[ram[!one.head, c(2,3), drop=FALSE]] <- P.pat[ram[!one.head, c(3,2), drop=FALSE]] <- 1
    C[observed, observed] <- R
    ref.indicators <- numeric(m)
    for (l in latent) {
        indicators <- A.pat[observed, l]
        if (any(indicators)){
            ref <- which.max(rowMeans(abs(R[observed[indicators], observed[indicators], drop=FALSE])))
            ref.indicators[l] <- observed[indicators][ref]
        }
        for (j in observed){
            C[j, l] <- C[l, j] <- if (!any(indicators)) runif(1, .3, .5)
            else {   
                numerator <- mean(abs(R[j, observed[indicators]]))
                denominator <- sqrt(mean(abs(R[observed[indicators], observed[indicators]])))
                c <- numerator/denominator
                sign(R[j, ref.indicators[l]])*abs(c)
            }
        }
    }
    for (l in latent){
        indicators.l <- A.pat[observed, l]
        for (k in latent){
            C[l, k] <- C[k, l] <-  if (l==k) 1 else { 
                #                 if (l == k){
                #                     if (!any(indicators.l)) runif(1, .3, .5)
                #                     else mean(abs(R[indicators.l, ref.indicators[l]]))
                #                 }
                #             else {
                indicators.k <- A.pat[observed, k]
                if ((!any(indicators.l)) | (!any(indicators.k))) runif(1, .3, .5) else {
                    numerator <- mean(abs(R[observed[indicators.l], observed[indicators.k]]))
                    denominator <- sqrt( mean(abs(R[observed[indicators.l], observed[indicators.l]]))
                                         * mean(abs(R[observed[indicators.k], observed[indicators.k]])))
                    c <- numerator/denominator
                    sign(R[ref.indicators[l], ref.indicators[k]])*abs(c)
                }
            }
        }
    }
    A <- matrix(0, m, m)
    for (j in 1:m){
        ind <- A.pat[j,]
        if (!any(ind)) next
        A[j, ind] <- solve(C[ind, ind]) %*% C[ind, j]
    }
    A[observed,] <- A[observed,]*matrix(s, n, m)
    A[,observed] <- A[,observed]/matrix(s, m, n, byrow=TRUE)
    C[observed,] <- C[observed,]*matrix(s, n, m)
    C[,observed] <- C[,observed]*matrix(s, m, n, byrow=TRUE)
    P <- (diag(m) - A) %*% C %*% t(diag(m) - A)
    P <- P.pat * P
    diag(P) <- abs(diag(P))
    for (par in 1:t){
        if (!is.na(start[par])) next
        posn <- par.posn[par]
        if (ram[posn, 1] == 1) start[par] <- A[ram[posn, 2], ram[posn, 3]]
        else start[par] <- P[ram[posn, 2], ram[posn, 3]]
        if (abs(start[par]) < tol) start[par] <- tol
    }
    if (debug){
        cat('\nStart values:\n')
        print(start)
        cat('\n')
    }
    start
}
