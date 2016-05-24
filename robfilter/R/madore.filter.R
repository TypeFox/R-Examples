#data(critvals)
#data(const)

madore.filter <- function(Y, byrow=FALSE, 
                             min.width=10, 
                             max.width=200,
                             test.sample.size=min.width/2,
                             width.search="geometric", 
                             rtr.size=min.width,
                             sign.level=0.1,
                             NA.sample.size=min.width, 
                             minNonNAs=min.width/2){

################
# Preparations #
################
    min.width <- round(min.width)
    max.width <- round(max.width)
    test.sample.size <- round(test.sample.size)
    rtr.size <- round(rtr.size)
    NA.sample.size <- round(NA.sample.size)
    minNonNAs <- round(minNonNAs)

##############################
# Stopping and Warning rules #
##############################
    if(missing(Y))
        stop("Input data set is missing with no default.\n")
        
    if(dim(Y)[1] < min.width)
        stop("Length of data set must be at least 'min.width' =", min.width, ".\n")
        
    if(!is.logical(byrow))
        stop("'byrow' must be either 'TRUE' or 'FALSE'.\n")
        
    if(!is.numeric(min.width))
        stop("'min.width' must be numeric.\n")
        
    if(!is.numeric(max.width))
        stop("'max.width' must be numeric.\n")
        
    if(max.width < min.width)
        stop("'max.width' must not be smaller than 'min.width'.\n")
        
    if(min.width < 10){
        min.width <- 10
        warning("'min.width' must be at least 10; 'min.width' is set to 10.\n")
        }
        
    if(!is.numeric(test.sample.size))
        stop("'test.sample.size' must be numeric.\n")
        
    if(test.sample.size < 5){
        test.sample.size <- 5
        warning("'test.sample.size' must be at least 5; 'test.sample.size' is set to 5.\n")
        }
             
    if(width.search!="linear" & width.search!="binary" & width.search!="geometric")
        stop("'width.search' must be a character, either 'linear', 'binary' or 'geometric'.\n")
        
    if(!is.numeric(rtr.size))
        stop("'rtr.size' must be numeric.\n")
        
    if(rtr.size < 0){
        rtr.size <- 0
        warning("'rtr.size' must be at least 0; 'rtr.size' is set to 0.\n")
        }
        
    if(!is.numeric(NA.sample.size))
        stop("'NA.sample.size' must be numeric.\n")
        
    if(NA.sample.size < 5){
        NA.sample.size <- 5
        warning("'NA.sample.size' must be at least 5; 'NA.sample.size' is set to 5.\n")
        }
        
    if(!is.numeric(minNonNAs))
        stop("'minNonNAs' must be numeric.\n")
        
    if(minNonNAs < 5){
        minNonNAs <- 5
        warning("'minNonNAs' must be at least 5; 'minNonNAs' is set to 5.\n")
    }

    if(sign.level <= 0 | sign.level > 0.5)
        stop("sign.level must be a value in (0,0.5)")

######################
# Internal functions #
######################

    get.test.sample.size <- function(n, test.sample.size){
        if(test.sample.size > n/2) return(floor(n/2))
        if(test.sample.size <= n/2) return(test.sample.size)
    }

    get.critval <- function(n, test.sample.size){
        if(n <= 600 & n > 10 & sign.level==0.1){
            return(as.numeric(critvals[n, test.sample.size]))  
        } else {
            return(2 * qhyper(1 - sign.level/2, floor(n/2), floor(n/2), test.sample.size) - test.sample.size)
        }
    }

    get.TS <- function(res){
        return(abs(sum(sign(res), na.rm = TRUE)))
    }

    OGK.qn <- function(X.mat, cqn){
        k <- dim(X.mat)[2]
        Qns <-  apply(X.mat, 2, Qn, cqn)
        too.small <- which(Qns < 0.02)
        Qns[too.small] <- 0.02
        D.mat <- diag(Qns)
        Y.mat <- X.mat %*% solve(D.mat)
        R.mat <- diag(1,k)
        for (i in 1:(k-1)){
            for (j in (i+1):k){
                R.mat[i,j] <- R.mat[j,i] <- 0.25*( Qn(X.mat[,i]+X.mat[,j],cqn)^2 - (Qn(X.mat[,i]-X.mat[,j],cqn)^2) )
            }
        }
        E.mat <- eigen(R.mat)$vectors
        A.mat <- D.mat %*% E.mat
        Z.mat <- X.mat %*% solve(t(A.mat))
        Gamma.diag <- apply(Z.mat, 2, Qn)
        too.small <- which(Gamma.diag < 0.02)
        Gamma.diag[too.small] <- 0.02
        Gamma.mat <- (diag(Gamma.diag))^2
        OGK <- A.mat %*% Gamma.mat %*% t(A.mat) 
        return(OGK)
    }

    LS.reg <- function(x, res, dat, cqn){
        di.qn  <- mahalanobis(res, center=rep(0,dim(dat)[2]), cov=OGK.qn(res, cqn), tol=2.220446e-16)
        dnqn   <- qchisq(0.95, dim(dat)[2]) * median(di.qn)/ qchisq(0.5, dim(dat)[2])
        kov    <- cov(cbind(x,dat)[di.qn <= dnqn,])
        beta   <- kov[1,2:(dim(dat)[2]+1)] / kov[1,1]
        sign.level  <- apply(dat[di.qn <= dnqn,], 2, mean) - (beta * mean(x))
        hat.y  <- sign.level + beta*(length(x))
        return(hat.y)
    }
    
    get.B.t <- function(y){
        n <- length(y)
        j <- 1:n
        B.t <- matrix(NA,n,n)
        for (i in j) {
            B.t[j[j<i],i] <- B.t[i,j[j<i]] <- (y[i] - y[j[j<i]]) / (i - j[j<i])
        }
        return(B.t)
    }
    
    get.beta.RM <- function(B.t){
        beta.i <- apply(B.t,1,median, na.rm=T)
        beta.RM <- median(beta.i, na.rm=T)
        return(beta.RM)
    }
    
    get.mu.RM <- function(y,beta.RM){
        n <- length(y)
        j <- 1:n
        mu.RM <- median(y - (beta.RM*(j-n)), na.rm=T)
        return(mu.RM)
    }

    aoRM.function <- function(x, B){
        n <- length(x)
        if(length(which(!is.na(x[(n-NA.sample.size+1):n]))) < minNonNAs){
        # there are not enough recent observations
            n.t <- NA
        } else { # there are enough recent observations
            beta.RM <- get.beta.RM(B)
            mu.RM <- get.mu.RM(x,beta.RM)
            x.hat <- mu.RM + beta.RM*((1:n)-n)
            res <- res <- x - x.hat
            test.sample.size2  <- get.test.sample.size(n, test.sample.size)
            TS <- get.TS(res[(n - test.sample.size2 + 1):n])
            if(TS > get.critval(n, test.sample.size2)){
                if(n > min.width){
                    direction <- -1
                } else {
                direction <- 0
                }
            } else {
                direction <- 0
            }
            if(direction==-1){
                if(width.search=="linear"){
                    n.t <- n
                    repeat{
                        if(n.t == min.width) break
                        n.t <- n.t-1
                        B.new <- B[(n-n.t+1):n, (n-n.t+1):n]
                        beta.RM <- get.beta.RM(B.new)
                        mu.RM <- get.mu.RM(x,beta.RM)
                        x.hat <- mu.RM + beta.RM*((1:n)-n)
                        res <- res <- x - x.hat
                        test.sample.size2  <- get.test.sample.size(n.t, test.sample.size)
                        TS <- get.TS(res[(n.t - test.sample.size2 + 1):n.t])
                        if(TS <= get.critval(n.t, test.sample.size2)) break
                    }
                }
                if(width.search=="binary"){
                    n.l <- min.width
                    n.u <- n
                    n.t <- n.l + floor((n.u-n.l)/2)
                }
                if(width.search=="geometric"){
                    j <- 1
                    repeat{
                        n.l <- n-(2^j)
                        if(n.l<=min.width){
                            n.l <- min.width
                            n.u <- n-2^(j-1)
                            n.t <- n.l + floor((n.u-n.l)/2)
                            break
                        }
                        B.new <- B[(n-n.l+1):n, (n-n.l+1):n]
                        x.win <- x[(n-n.l+1):n]
                        
                        beta.RM <- get.beta.RM(B.new)
                        mu.RM <- get.mu.RM(x.win,beta.RM)
                        x.hat <- mu.RM + beta.RM*((1:n.l)-n.l)
                        res <- x.win - x.hat
                        test.sample.size2  <- get.test.sample.size(n.l, test.sample.size)
                        TS <- get.TS(res[(n.l - test.sample.size2 + 1):n.l])
                        if(TS > get.critval(n.l, test.sample.size2)){
                            j <- j+1
                        } else {
                            n.u <- n-2^(j-1)
                            n.t <- n.l + floor((n.u-n.l)/2)
                            break
                        }
                    }
                }
                if(width.search!="linear"){
                    repeat{
                        x.win <- x[(n-n.t+1):n]
                        B.new <- B[(n-n.t+1):n, (n-n.t+1):n]
                        beta.RM <- get.beta.RM(B.new)
                        mu.RM <- get.mu.RM(x.win,beta.RM)
                        x.hat <- mu.RM + beta.RM*((1:n.t)-n.t)
                        res <- x.win - x.hat
                        test.sample.size2  <- get.test.sample.size(n.t, test.sample.size)
                        TS <- get.TS(res[(n.t - test.sample.size2 + 1):n.t])
                        if(n.t==n.l){
                            break
                        }
                        if(n.t==n.u & TS <= get.critval(n.t, test.sample.size2)){
                            break
                        }
                        if(n.t==n.u & TS > get.critval(n.t, test.sample.size2)){
                            n.t <- n.l
                            break
                        }
                        if(TS > get.critval(n.t, test.sample.size2)){
                            n.u <- n.t
                            n.t <- n.l + floor((n.u-n.l)/2)
                        } else {
                            n.l <- n.t
                            n.t <- n.l + ceiling((n.u-n.l)/2)
                        }
                    }
                }
            }
            if(direction==0){
                n.t <- n
            }
        }
        return(n.t)
    }

####################
# Internal objects #
####################
    if(byrow==TRUE) Y <- t(Y)
    N <- dim(Y)[1]
    K <- dim(Y)[2]
    y <- matrix(NA,N,K)
    for(k in 1:K){
        y[,k] <- as.vector(Y[,k])
    }
    Y <- y
    if(!is.numeric(Y))
        stop("Data set must be numeric.\n")
    n <- round(min.width)
    signals <- widths <- matrix(NA, nrow=N, ncol=K)
    ov.width <- rep(NA,N)
    B <- NULL
    for(k in 1:K){
        B[[k]] <- get.B.t(Y[1:n,k])
    }
    indiv.widths <- rep(NA,K)

###################
# Main  Algorithm #
###################
    for(i in n:N){
        # apply aoRM to obtain individual window widths
        Y.win <- as.matrix(Y[(i-n+1):i,])
        for(k in 1:K){
            indiv.widths[k] <- aoRM.function(x=Y.win[,k], B=B[[k]])
        }
        widths[i,] <- indiv.widths
        if(any(!is.na(indiv.widths))){ # there is at least one variable that offers enough recent observations
            ov.width[i] <- n.t <- min(indiv.widths,na.rm=TRUE)
            position <- which(!is.na(indiv.widths))
            Y.win2 <- Y[(i-n.t+1):i,position]
            K.t <- length(position)
            if(K.t == 1){
                index <- (n-indiv.widths[position]+1):n
                beta.RM <- get.beta.RM(B[[position]][index,index])
                mu.RM <- get.mu.RM(Y.win2, beta.RM)
                signals[i,position] <- mu.RM
            }
            if(K.t > 1){
                res <- matrix(NA, ncol=K.t, nrow = ov.width[i])
                for(k in 1:K.t){
                    nk <- dim(B[[position[k]]])[1]
                    beta.RM <- get.beta.RM(B[[position[k]]][(nk-n.t+1):nk,(nk-n.t+1):nk])
                    mu.RM <- get.mu.RM(Y.win2[,k],beta.RM)
                    x.hat <- mu.RM + beta.RM*((1:n.t)-n.t) 
                    j <- which(is.na(Y.win2[,k]))
                    if(length(j) > 0){
                        Y.win2[j,k] <- x.hat[j]
                    }
                    res[,k] <- Y.win2[,k] - x.hat
                }
                x <- 1:ov.width[i]
                cqn <- const$const[which.min(abs(const$n - ov.width[i]))]
                signals[i,position] <- LS.reg(x, res, Y.win2, cqn)
                if(rtr.size>0){
                    if(rtr.size > ov.width[i]){
                        for(k in 1:K.t){
                            if(signals[i,position[k]] > max(Y.win[,position[k]], na.rm = TRUE)){
                                signals[i,position[k]] <- max(Y.win[,position[k]], na.rm = TRUE)
                            }
                            if(signals[i,position[k]] < min(Y.win[,position[k]], na.rm = TRUE)){
                                signals[i,position[k]] <- min(Y.win[,position[k]], na.rm = TRUE)
                            }
                        }
                    } else {
                        for(k in 1:K.t){
                            if(all(is.na(Y.win[(n-rtr.size+1):n,k]))){
                                if(signals[i,position[k]] > max(Y.win[,position[k]], na.rm = TRUE)){
                                    signals[i,position[k]] <- max(Y.win[,position[k]], na.rm = TRUE)
                                }
                                if(signals[i,position[k]] < min(Y.win[,position[k]], na.rm = TRUE)){
                                    signals[i,position[k]] <- min(Y.win[,position[k]], na.rm = TRUE)
                                }
                            } else {
                                if(signals[i,position[k]] > max(Y.win[(n-rtr.size+1):n,position[k]], na.rm = TRUE)){
                                    signals[i,position[k]] <- max(Y.win[(n-rtr.size+1):n,position[k]], na.rm = TRUE)
                                }
                                if(signals[i,position[k]] < min(Y.win[(n-rtr.size+1):n,position[k]], na.rm = TRUE)){
                                    signals[i,position[k]] <- min(Y.win[(n-rtr.size+1):n,position[k]], na.rm = TRUE)
                                }
                            }
                        }
                    }
                }
            }
            # update-step
            if(i != N){
                if(ov.width[i] < max.width){
                    n <- ov.width[i]+1
                }
                for(k in 1:K){
                    nk <- dim(B[[k]])[1]
                    b.new  <- (Y[i+1,k] - Y[(i-n+2):i,k]) / (i+1 - (i-n+2):i)
                    B[[k]] <- cbind(rbind(B[[k]][(nk-n+2):nk, (nk-n+2):nk], b.new), c(b.new,NA))
                }
            }    
        } else { # there is no variable that offers enough recent observations
            # update-step
            if(i != N){     
                n <- min.width       
                for(k in 1:K){
                    nk <- dim(B[[k]])[1]
                    b.new  <- (Y[i+1,k] - Y[(i-n+2):i,k]) / (i+1 - (i-n+2):i)
                    B[[k]] <- cbind(rbind(B[[k]][(nk-n+2):nk, (nk-n+2):nk], b.new), c(b.new,NA))   
                }  
            }
        } 
    }
    result <- list(signals=signals, widths=widths, overall.width=ov.width, 
                   Y=Y, byrow=byrow, min.width=min.width,
                   max.width=max.width, test.sample.size=test.sample.size, width.search=width.search,
                   rtr.size=rtr.size, NA.sample.size=NA.sample.size, minNonNAs=minNonNAs)
    return(structure(result, class="madore.filter"))
}

##################
# Default output #
##################
print.madore.filter <- function(x, ...){
    N <- dim(x$signals)[1]
    cat('$signals \n')
    if(N <= 100){
        print(x$signals, ...)
    } else {
        print(x$signals[1:(x$min.width + 9),])
        cat('Only the first 10 signal estimations are printed.\n')
    }
}

################
# Default plot #
################
plot.madore.filter <- function(x, ...){
    N <- dim(x$Y)[1]
    k <- dim(x$Y)[2]
    plot(rep(NA,N), main="Multivariate time series and signal extraction by madore.filter",
    type='l', xlab='time', ylab='', ylim=c(min(x$Y, na.rm=TRUE),max(x$Y, na.rm=TRUE)), las=1,...)
    for(i in 1:k){
        lines(x$Y[,i])
        lines(x$signals[,i], col=2, lwd=1)
    }
    legend(x="topleft", bty="n", legend=c("Data", "Signal extraction"), lty=c(1,1), col=c(1,2), lwd=c(1,1))
}
