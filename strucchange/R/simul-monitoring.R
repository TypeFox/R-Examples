monitorMECritvalData <- function(n, end=10, frequency=1000, h=1,
                           verbose=FALSE){
    logPlus <- function(x) ifelse(x<=exp(1),1,log(x))
    
    dn <- list(NULL, c("max","range"), as.character(h), as.character(end))
    
    z <- array(0, dim=c(n, 2, length(h), length(end)), dimnames=dn)
    e <- max(end)
    t <- seq(1+1/frequency,e,by=1/frequency)
    
    l <- logPlus(t)
    
    H <- h*frequency
    
    for(k in 1:n){
        if(verbose){
            cat("  Loop:    ",k,"\r")
        }
        b <- as.vector(e1071::rbridge(e,frequency))
        for(hh in 1:length(h)){
            bh <- c(rep(0, length=H[hh]), diff(b, lag=H[hh]))
            bh <- bh[(frequency+1):length(bh)]
            bh <- bh/sqrt(2*l)
            
            for(ee in 1:length(end)){
                index <- 1:((end[ee]-1)*frequency)
                z[k,1,hh,ee] <- max(abs( bh[index] ))
                z[k,2,hh,ee] <- max(bh[index]) - min(bh[index])
            }
        }
    }
    cat("                               \r")
    z
}

monitorMECritval <- function(x, probs=c(0.9,0.95))
{    
    dx <- dim(x)
    lp <- length(probs)

    z<-array(0,c(dx[3],dx[4],lp,2))
    dimnames(z) <- list(dimnames(x)[[3]],dimnames(x)[[4]],
                        probs, c("max", "range"))
    for(k in 1:dx[3]){
        for(l in 1:dx[4]){
            for(m in 1:lp){
                z[k,l,,1] <- quantile(x[,"max",k,l],probs)
                z[k,l,,2] <- quantile(x[,"range",k,l],probs)
            }
        }
    }
    
    z
}

monitorRECritvalData <- function(n, end=10, frequency=1000,
                                 verbose=FALSE, linear=FALSE){
    
    dn <- list(NULL, as.character(end))
    
    z <- array(0, dim=c(n, length(end)), dimnames=dn)
    e <- max(end)
    t <- seq(1+1/frequency,e,by=1/frequency)
    
    for(k in 1:n){
        if(verbose){
            cat("  Loop:    ",k,"\r")
        }
        b <- as.vector(e1071::rbridge(e,frequency))[-(1:frequency)]
        if(linear)
#            b <- b/pmax(1, 0.5*t+0.25)
            b <- b / t
        else
            b <- b^2/(t*(t-1)) - log(t/(t-1))
        
        for(ee in 1:length(end)){
            index <- 1:((end[ee]-1)*frequency)
            z[k,ee] <- max(abs(b[index]))
        }
    }
    cat("                               \r")
    if(linear)
        attr(z, "square") <- FALSE
    else
        attr(z, "square") <- TRUE
    z
}

monitorRECritval <- function(x, probs=c(0.9,0.95))
{    
    dx <- dim(x)
    lp <- length(probs)

    z<-array(0, dim=c(dx[2],lp))
    dimnames(z) <- list(dimnames(x)[[2]], probs)

    for(k in 1:dx[2]){
            z[k,] <- quantile(x[,k],probs)
    }

    if(attr(x, "square"))
        z <- sqrt(z)

    z
}

