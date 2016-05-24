#data(dfs)
#data(var.n)

scarm.filter <- function(time.series, 
                         right.width=30, 
                         min.left.width=right.width, 
                         min.width=floor(right.width/3), 
                         max.width=200, 
                         sign.level=0.001,
                         bound.noise.sd=0.01,
                         rtr=TRUE,
                         autocorrelations="automatic"
                         ){

################
# Preparations #
################
    r <- round(right.width)
    ell.min <- round(min.left.width)
    n.min <- round(min.width)
    n.max <- round(max.width)
    length.series <- length(time.series) 

##############################
# Stopping and Warning rules #
##############################

    if(missing(time.series))
        stop("The input data vector is missing with no default.\n")
        
    if(!is.numeric(time.series))
        stop("The data vector must be numeric.\n")
        
    if(length.series < min.width)
        stop("Length of time series must be at least 'min.width' =", min.width, ".\n")

    if(min.width < 5){
        min.width <- 5
        warning("'min.width' must be an integer >=5; 'min.width' is set to 5.\n")
    }
    
    if(right.width < 5){
        right.width <- 5
        warning("'right.width' must be an integer >=5; 'right.width' is set to 5.\n")
    }    
        
    if(min.left.width < 5){
        min.left.width <- 5
        warning("'min.left.width' must be an integer >=5; 'min.left.width' is set to 5.\n")
    }  
    
    if(right.width > min.left.width){
        min.left.width <- right.width
        warning("'min.left.width' must not be smaller than 'right.width'; 'min.left.width' is set to 'right.width' = ", right.width,".\n")
    }

    if(sign.level <= 0 | sign.level > 0.5)
        stop("'sign.level' must be a value in (0,0.5)")
        
    if(max.width < min.width)
        stop("'max.width' must not be smaller than 'min.width'.\n")
        
    if(max.width < right.width + min.left.width | max.width%%1!=0)
        stop("'max.width' must not be smaller than 'right.width'+'min.left.width' = ", right.width+min.left.width,".\n")
    
    if(bound.noise.sd <= 0){
        bound.noise.sd <- 0.01
        warning("'bound.noise.sd' must be a value >0; 'bound.noise.sd' is set to 0.01.\n")
    }
    
    if(!is.logical(rtr)){
        stop("'rtr' must be either TRUE or FALSE.\n")
    }
    
    if(all(autocorrelations!=c("high.positive", "moderate.positive", "small.positive", "no", "small.negative", "moderate.negative", "high.negative", "automatic"))){
        stop("'autocorrelations' must be either 'high.positive', 'moderate.positive', 'small.positive', 'no', 'small.negative', 'moderate.negative', 'high.negative' or 'automatic'.\n")
    }
                                          
#############################################
# Internal functions and required data sets #
#############################################

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
    
    get.B1 <- function(y, ell){
        n <- length(y)
        r <- n - ell
        I.left <- 1:ell
        I.right <- (ell+1):n
        B1 <- matrix(NA,ell,r)
        for (i in I.left) {
            B1[i,I.right-ell] <- (y[i] - y[I.right]) / (i - I.right)
        }
        return(B1)
    }

    # Q scale estimator with bias-correction c.q depending on chosen value for 'autocorrelations'
    if(autocorrelations=="high.negative"){
        c.q.sim   <- const.Q[,1]
        c.q.model <- function(n){
            c.q.sim[300]
        }
    }
    if(autocorrelations=="moderate.negative"){
        c.q.sim   <- const.Q[,2]
        c.q.model <- function(n){
             c.q.sim[300]
        }
    }
    if(autocorrelations=="small.negative"){
        c.q.sim   <- const.Q[,3]
        c.q.model <- function(n){
             c.q.sim[300]
        }
    }
    if(autocorrelations=="no"){
        c.q.sim   <- const.Q[,4]
        c.q.model <- function(n){
            1.21 * (n/(n+0.44))
        }
    }
    if(autocorrelations=="small.positive"){
        c.q.sim   <- const.Q[,5]
        c.q.model <- function(n){
             c.q.sim[300]
        }
    }
    if(autocorrelations=="moderate.positive"){
        c.q.sim   <- const.Q[,6]
        c.q.model <- function(n){
             c.q.sim[300]
        }
    }
    if(autocorrelations=="high.positive"){
        c.q.sim   <- const.Q[,7]
        c.q.model <- function(n){
             c.q.sim[300]
        }
    }

    get.Q.adj <- function(y, phi.est){
        n <- length(y)
        m <- length(na.omit(y))
        if(m < 3)
            stop("data vector must contain at least 3 observations")
        x <- which(!is.na(y))
        y <- y[x]
        if(autocorrelations=="automatic"){
            autocor.type <- which.min(abs(phi.est-c(-0.9,-0.6,-0.3,0,0.3,0.6,0.9)))
            c.q.sim   <- const.Q[,autocor.type]
            c.q.model <- function(n){
                c.q.sim[300]
            }
        }
        if(m > 4){
            if(m <= 150){
                c.q <- c.q.sim[m]
            } else {
                c.q <- c.q.model(m)
            }
        } else {
            c.q <- 1
        }
        h <- numeric(m-2)
        for(j in 1:(m-2)){
            h[j] <- abs(y[j+1] - y[j] - (x[j+1] - x[j])*((y[j+2]-y[j])/(x[j+2]-x[j])))
        } 
        h <- sort(h)
        out <- c.q * sort(h)[max(1,floor(0.5*(m-2)))]
        return(max(out,bound.noise.sd))
    }

    # degrees of freedom for t-distribution to obtain critical values for the SCARM test statistic
    get.critval <- function(ell, r, sign.level){
        if(ell <= 100){
            ell <- ell - ell%%5
            r <- r - r%%5
            degree.of.freedom <- dfs[ell/5,r/5]
            critval <- qt(p=1-(sign.level/2), df=degree.of.freedom)
        } else {
            critval <- qnorm(p=1-(sign.level/2))
            
        }
        return(critval)
    }
    
    
    if(autocorrelations=="high.negative"){
        emp.var <- var.n[,1]
        v.model <- function(n){
            out <- 0.0003514981 + (-0.0629949278/n) + (1.6446684583/(n^2)) + (-3.1859306183/(n^3))
            return(out)
        }
    }
    if(autocorrelations=="moderate.negative"){
        emp.var <- var.n[,2]
        v.model <- function(n){
            out <- 0.0002560433 + (-0.0434806450/n) + (0.9224130380/(n^2)) + (5.4903584072/(n^3))
            return(out)
        }
    }     
    if(autocorrelations=="small.negative"){
        emp.var <- var.n[,3]
        v.model <- function(n){
            out <- 0.0002385982 + (-0.0395730416/n) + (0.7798496595/(n^2)) + (9.7595002645/(n^3))
            return(out)
        } 
    } 
    if(autocorrelations=="no"){
        emp.var <- var.n[,4]
        v.model <- function(n){
            out <- 4.77e-07 + (17.71*(1/n^3))
            return(out)
        }
    }
    if(autocorrelations=="small.positive"){
        emp.var <- var.n[,5]
        v.model <- function(n){
            out <- 0.000473727 + (-0.086723762/n) + (2.442400804/(n^2)) + (10.578915504/(n^3))
            return(out)
        } 
    }        
    if(autocorrelations=="moderate.positive"){
        emp.var <- var.n[,6]
        v.model <- function(n){
            out <- 0.0007355669 + (-0.1485876684/n) + (5.5086877566/(n^2)) + (-5.3249305456/(n^3))
            return(out)
        }
    }   
    if(autocorrelations=="high.positive"){
        emp.var <- var.n[,7]
        v.model <- function(n){
            out <- -0.0003318919 + (0.0476670483/n) + (2.2422001348/(n^2)) + (-5.8731008014/(n^3))
            return(out)
        }
    }

    # this function delivers the scarm test statistic
    scarm.test.statistic <- function(y, B, ell, r){
        n <- ell+r
        left.sample <- y[1:ell]
        right.sample <- y[(ell+1):n]
        B.t.left <- B[1:ell,1:ell]
        B.t.right <- B[(ell+1):n,(ell+1):n]
        beta.left <- get.beta.RM(B.t.left)
        beta.right <- get.beta.RM(B.t.right)
        d.t <- beta.left - beta.right
        if(autocorrelations=="automatic"){
            autocor.type <- which.min(abs(phi.est-c(-0.9,-0.6,-0.3,0,0.3,0.6,0.9)))
            emp.var <- var.n[,autocor.type]
            v.model <- function(n){
                emp.var[300]
            }
        }
        if(ell<=300){
            v.left <- emp.var[ell]
        } else {
            v.left <- v.model(ell)
        }
        if(r<=300){
            v.right <- emp.var[r]
        } else {
            v.right <- v.model(r)
        }
        var.t <- noise.sd.est^2 * (v.left+v.right)
        test.statistic <- d.t/sqrt(var.t)
        return(list(test.statistic=test.statistic, beta.left=beta.left, beta.right=beta.right))
    }

####################
# Internal objects #
####################
    signal.est <- slope.est <- slope.diff <- adapted.width <- scarm.statistic <- critvals <- noise.sd <- acf.lag.one <- rep(NA,length.series)  
    n <- n.min
    I.win <- 1:n
    y <- time.series[I.win]
    B <- get.B.t(y)

###################
# Main  Algorithm #
###################
    for(i in n.min:length.series){
        phi.est.sample.data <- time.series[(max(i-n.max+1,1)):i]
        phi.est.sample.signal <- signal.est[(max(i-n.max+1,1)):i]
        phi.est.sample <- na.omit(phi.est.sample.data-phi.est.sample.signal)
        if(length(phi.est.sample)>=n.min){
            phi.est <- acf(phi.est.sample, plot=F)$acf[2]
            if(is.na(phi.est)){
                phi.est <- 0
            }
        } else {
            phi.est <- 0
        }
        if(all(is.na(y[(n-n.min+1):n])) || length(which(!is.na(y[(n-(min(r,n))+1):n]))) < n.min){
            # there are not enough recent observations!
            mu.est <- NA
            beta.est <- NA
            noise.sd.est <- NA
            # Update step
            m <- min(n,ell.min+r)
            adapted.width[i] <- m
            if(i != length.series){
                b.new <- (time.series[i+1] - time.series[i+c((-m + 1):0)]) / (i+1 - (i+c((-m + 1):0)))
                B <- B[(n-m+1):n,(n-m+1):n]     
                B     <- cbind(rbind(B, b.new), c(b.new, NA))
                y     <- time.series[(i - m + 1):(i + 1)]
                n     <- m + 1
            }
        } else {
            if(n >= ell.min+r){
            # is n large enough for test application?
                ell <- n - r
                if(length(which(!is.na(y[1:ell])))>=ell.min/2 && length(which(!is.na(y[(ell+1):n])))>=r/2){
                # are there enough observations for testing? 
                    noise.sd.est <- get.Q.adj(y, phi.est=phi.est)
                    critval <- get.critval(ell=ell, r=r, sign.level=sign.level)
                    scarm.test <- scarm.test.statistic(y=y, B=B, ell=ell, r=r)
                    scarm.statistic[i] <- scarm.test$test.statistic
                    slope.diff[i] <- scarm.test$beta.left - scarm.test$beta.right
                    critvals[i] <- critval
                    noise.sd[i] <- noise.sd.est
                    if(abs(scarm.test$test.statistic) > critval){
                    # does test reject the null hypothesis?
                        B <- B[(n-n.min+1):n,(n-n.min+1):n]            
                        beta.est <- get.beta.RM(B)
                        mu.est <- get.mu.RM(beta.RM=beta.est, y=y[(n-n.min+1):n])
                        n <- n.min
                        if(rtr==TRUE){
                            rtr.sample <- time.series[(i-n.min+1):i]
                            min.rtr.sample <- min(rtr.sample, na.rm=T)
                            max.rtr.sample <- max(rtr.sample, na.rm=T)
                            mu.est <- max(min.rtr.sample, min(mu.est, max.rtr.sample))
                        }
                        if(i != length.series) {
                            b.new <- (time.series[i+1] - time.series[i+c((-n + 1):0)]) / (i+1 - (i+c((-n + 1):0)))
                            B <- cbind(rbind(B, b.new), c(b.new,NA))
                            adapted.width[i] <- n
                            y <- time.series[(i - n + 1):(i + 1)]
                            n <- n+1
                        }
                    } else {
                    # difference of RM slopes is NOT too large => do not decrease time window
                        beta.est <- get.beta.RM(B)
                        mu.est <- get.mu.RM(beta.RM=beta.est, y=y)                    
                        adapted.width[i] <- n
                        if(i != length.series){
                            if(n == n.max){ 
                            # is window width maximal?
                                b.new     <- (time.series[i+1] - time.series[i+c((-n + 2):0)]) / (i+1 - (i+c((-n + 2):0)))
                                B     <- cbind(rbind(B[-1, -1], b.new), c(b.new,NA))
                                y       <- time.series[(i - n + 2):(i + 1)]
                            } else {
                                b.new     <- (time.series[i+1] - time.series[i+c((-n + 1):0)]) / (i+1 - (i+c((-n + 1):0)))
                                B       <- cbind(rbind(B, b.new), c(b.new, NA))
                                y       <- time.series[(i - n + 1):(i + 1)]
                                n       <- n + 1
                            }
                        }
                        if(rtr==TRUE){
                            rtr.sample <- time.series[(i-n.min+1):i]
                            min.rtr.sample <- min(rtr.sample, na.rm=T)
                            max.rtr.sample <- max(rtr.sample, na.rm=T)
                            mu.est <- max(min.rtr.sample, min(mu.est, max.rtr.sample))
                        }    
                    }
                } else {
                # there are not enough observations for testing!
                    B <- B[(n-ell.min-r+1):n,(n-ell.min-r+1):n]
                    beta.est <- get.beta.RM(B)
                    mu.est <- get.mu.RM(beta.RM=beta.est, y=y[(n-ell.min-r):n])
                    #regression.line <- mu.est + beta.est*((1-r):0)
                    #res <- y[(n-r+1):n] - regression.line
                    noise.sd.est <- get.Q.adj(y, phi.est=phi.est)
                    n <- ell.min+r
                    if(rtr==TRUE){
                        rtr.sample <- time.series[(i-n.min+1):i]
                        min.rtr.sample <- min(rtr.sample, na.rm=T)
                        max.rtr.sample <- max(rtr.sample, na.rm=T)
                        mu.est <- max(min.rtr.sample, min(mu.est, max.rtr.sample))
                    }
                    if(i != length.series) {
                        b.new <- (time.series[i+1] - time.series[i+c((-n + 1):0)]) / (i+1 - (i+c((-n + 1):0)))
                        B <- cbind(rbind(B, b.new), c(b.new,NA))
                        adapted.width[i] <- n
                        y <- time.series[(i - n + 1):(i + 1)]
                        n <- n+1
                    }  
               }   
           } else {
           # n is not large enough for test application!
                beta.est <- get.beta.RM(B)
                mu.est <- get.mu.RM(y=y, beta.RM=beta.est)
                #regression.line <- mu.est + beta.est*((1-n):0)
                #res <- y - regression.line
                noise.sd.est <- get.Q.adj(y, phi.est=phi.est)
                adapted.width[i] <- n
                if(rtr==TRUE){
                    rtr.sample <- time.series[(i-n.min+1):i]
                    min.rtr.sample <- min(rtr.sample, na.rm=T)
                    max.rtr.sample <- max(rtr.sample, na.rm=T)
                    mu.est <- max(min.rtr.sample, min(mu.est, max.rtr.sample))
                }
                if(i != length.series) {
                    if(n == n.max){ 
                    # is window width maximal?
                        b.new     <- (time.series[i+1] - time.series[i+c((-n + 2):0)]) / (i+1 - (i+c((-n + 2):0)))
                        B     <- cbind(rbind(B[-1, -1], b.new), c(b.new,NA))
                        y       <- time.series[(i - n + 2):(i + 1)]
                    } else {
                        b.new     <- (time.series[i+1] - time.series[i+c((-n + 1):0)]) / (i+1 - (i+c((-n + 1):0)))
                        B       <- cbind(rbind(B, b.new), c(b.new, NA))
                        y       <- time.series[(i - n + 1):(i + 1)]
                        n       <- n + 1
                    }
                }       
            }
        }
        signal.est[i] <- mu.est
        slope.est[i] <- beta.est
        noise.sd[i] <- noise.sd.est
        acf.lag.one[i] <- phi.est
    }
    result <- list(signal.est=signal.est, slope.est=slope.est, adapted.width=adapted.width, 
                   test.statistic=scarm.statistic, critvals=critvals, 
                   noise.sd=noise.sd, slope.diff=slope.diff, acf.lag.one=acf.lag.one,
                   time.series=time.series, right.width=r, min.left.width=ell.min, min.width=n.min, max.width=n.max, 
                   sign.level=sign.level, bound.noise.sd=bound.noise.sd, 
                   rtr=rtr, autocorrelations=autocorrelations)
    return(structure(result, class = "scarm.filter"))
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Default output
print.scarm.filter <- function(x, ...) {
    length.series <- length(x$signal)
    cat('$signal \n')  
      if(length.series <= 100){
      print(x$signal,...)
    } else {
      print(x$signal[1:50])
      cat('Only the first 50 signal estimations are printed.\n')    
    }
}

## # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
## Default plot
plot.scarm.filter <- function(x, info=T,...){
    if(info==T){
        par(mfrow=c(2,1))
        plot(x$time.series, pch=20, cex=0.5, main="Time series and SCARM signal extraction", xlab="time", ylab="", ...)
        lines(x$time.series)
        lines(x$signal, col=2, lwd=1)
        legend(x="topleft", legend=c("Data", "SCARM signal extraction"),lty=c(1, 1),col=c("black", "red"),lwd=c(1, 1), bty="n")
        plot(x$adapted.width, type="l",  main="Adapted window widths", xlab="time", ylab="", ...)
    } else {
        par(mfrow=c(1,1))        
        plot(x$time.series, pch=20, cex=0.5, main="Time series and SCARM signal extraction", xlab="time", ylab="", ...)
        lines(x$time.series)
        lines(x$signal, col=2, lwd=1)
        legend(x="topleft", legend=c("Data", "SCARM signal extraction"),lty=c(1, 1),col=c("black", "red"),lwd=c(1, 1), bty="n")
    }
}
