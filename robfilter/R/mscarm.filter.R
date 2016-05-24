#data(dfs)
#data(const.Q)
#data(var.n)

mscarm.filter <- function(time.series,
                          # input arguments of the univariate SCARM
                          right.width=30,
                          min.left.width=right.width,
                          min.width=floor(right.width/3),
                          max.width=200,
                          sign.level=0.001,
                          bound.noise.sd=0.01,
                          rtr=TRUE,                          
                          autocorrelations="automatic",
                          # input arguments of the SSM
                          c.bound=3,
                          r.bound=0
                          ){

    #######################################
    # Stopping and rules and preparations #
    #######################################
    if(missing(time.series))
        stop("The input data is missing with no default.\n")
    if(!is.numeric(time.series))
        stop("Data must be numeric.\n")
    if(!is.matrix(time.series))
        stop("'time.series' must be a matrix.\n")    
    length.series <- dim(time.series)[1]
    dim.series <- dim(time.series)[2]
    if(dim.series<2)
        stop("'time.series' must be a multivariate time series of dimension >=2.\n")
    for(i in 1:dim.series){
        if(all(is.na(time.series[,i])))
            stop("At least one variable has only missing values.\n")
    }
    if(length.series < dim.series){
        time.series <- t(time.series)
        length.series <- dim(time.series)[1]
        dim.series <- dim(time.series)[2]
    }
    if(is.matrix(time.series)==FALSE){
        time.series <- as.matrix(time.series)
    }
    r <- round(right.width)
    ell.min <- round(min.left.width)
    n.min <- round(min.width)
    n.max <- round(max.width)
    if(length.series < n.min)
        stop("Length of time series must be at least 'n.min' =", n.min, ".\n")
    if(n.min < 5){
        n.min <- 5
        warning("'n.min' must be an integer >=5; 'n.min' is set to 5.\n")
    }
    if(r < 5){
        r <- 5
        warning("'r' must be an integer >=5; 'r' is set to 5.\n")
    }  
    if(ell.min < 5){
        ell.min <- 5
        warning("'ell.min' must be an integer >=5; 'ell.min' is set to 5.\n")
    }  
    if(r > ell.min){
        ell.min <- r
        warning("'ell.min' must not be smaller than 'r'; 'ell.min' is set to 'r' = ", r,".\n")
    }
    if(sign.level <= 0 | sign.level > 0.5)
        stop("'sign.level' must be a value in (0,0.5)")
    if(n.max < n.min)
        stop("'n.max' must not be smaller than 'n.min'.\n")
    if(n.max < r + ell.min | n.max%%1!=0)
        stop("'n.max' must not be smaller than 'r'+'ell.min' = ", r+ell.min,".\n")
    if(bound.noise.sd <= 0){
        bound.noise.sd <- 0.01
        warning("'bound.noise.sd' must be a value >0; 'bound.noise.sd' is set to 0.01.\n")
    }
    if(c.bound <= 0){
        c.bound <- 3
        warning("'c.bound' must be a value >0; 'c.bound' is set to 3.\n")
    }
    if(r.bound < 0 | r.bound > 1){
        r.bound <- 0
        warning("'r.bound' must be a value in [0,1]; 'r.bound' is set to 0.\n")
    }
    if(!is.logical(rtr))
        stop("'rtr' must be either TRUE or FALSE.\n")
    if(all(autocorrelations!=c("high.positive", "moderate.positive", "small.positive", "no", "small.negative", "moderate.negative", "high.negative", "automatic")))
        stop("'autocorrelations' must be either 'high.positive', 'moderate.positive', 'small.positive', 'no', 'small.negative', 'moderate.negative', 'high.negative' or 'automatic'.\n")

    ######################
    # required functions #
    ######################
    
    # Q scale estimator with bias-correction c.q depending on chosen value for 'autocorrelations'
    if(autocorrelations=="high.negative"){
        c.q.sim   <- const.Q[,1]
        c.q.model <- function(n){
            c.q.sim[300]
        }
        emp.var <- var.n[,1]
        v.model <- function(n){
            out <- 0.0003514981 + (-0.0629949278/n) + (1.6446684583/(n^2)) + (-3.1859306183/(n^3))
            return(out)
        }
    }
    if(autocorrelations=="moderate.negative"){
        c.q.sim   <- const.Q[,2]
        c.q.model <- function(n){
             c.q.sim[300]
        }
        emp.var <- var.n[,2]
        v.model <- function(n){
            out <- 0.0002560433 + (-0.0434806450/n) + (0.9224130380/(n^2)) + (5.4903584072/(n^3))
            return(out)
        }
    }
    if(autocorrelations=="small.negative"){
        c.q.sim   <- const.Q[,3]
        c.q.model <- function(n){
             c.q.sim[300]
        }
        emp.var <- var.n[,3]
        v.model <- function(n){
            out <- 0.0002385982 + (-0.0395730416/n) + (0.7798496595/(n^2)) + (9.7595002645/(n^3))
            return(out)
        } 
    }
    if(autocorrelations=="no"){
        c.q.sim   <- const.Q[,4]
        c.q.model <- function(n){
            1.21 * (n/(n+0.44))
        }
        emp.var <- var.n[,4]
        v.model <- function(n){
            out <- 4.77e-07 + (17.71*(1/n^3))
            return(out)
        }
    }
    if(autocorrelations=="small.positive"){
        c.q.sim   <- const.Q[,5]
        c.q.model <- function(n){
             c.q.sim[300]
        }
        emp.var <- var.n[,5]
        v.model <- function(n){
            out <- 0.000473727 + (-0.086723762/n) + (2.442400804/(n^2)) + (10.578915504/(n^3))
            return(out)
        } 
    }
    if(autocorrelations=="moderate.positive"){
        c.q.sim   <- const.Q[,6]
        c.q.model <- function(n){
             c.q.sim[300]
        }
        emp.var <- var.n[,6]
        v.model <- function(n){
            out <- 0.0007355669 + (-0.1485876684/n) + (5.5086877566/(n^2)) + (-5.3249305456/(n^3))
            return(out)
        }
    }
    if(autocorrelations=="high.positive"){
        c.q.sim   <- const.Q[,7]
        c.q.model <- function(n){
             c.q.sim[300]
        }
        emp.var <- var.n[,7]
        v.model <- function(n){
            out <- -0.0003318919 + (0.0476670483/n) + (2.2422001348/(n^2)) + (-5.8731008014/(n^3))
            return(out)
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
    
    # covariance estimator OGK based on scale estimator Q.adj
    get.Ogk <- function(x, noise.sd.est.s){
       D.Ogk <- diag(noise.sd.est.s)
       Y.Ogk <- x%*%solve(D.Ogk)
       R.Ogk <- diag(1, nrow=ncol(x))
       Rij.Ogk <- function(x){
           ((get.Q.adj(apply(x,1,sum), phi.est=0))^2 -  (get.Q.adj(apply(x,1,diff), phi.est=0))^2)/4
       }
       comparisons <- t(combn(1:ncol(x),2))
       for(i in 1:nrow(comparisons)){
           column1 <- comparisons[i,1]
           column2 <- comparisons[i,2]
           R.Ogk[column1,column2] <-  R.Ogk[column2,column1] <- Rij.Ogk(cbind(Y.Ogk[,column1], Y.Ogk[,column2]))
       }
       Eigen.R.Ogk <- eigen(R.Ogk)
       E.Ogk <-  Eigen.R.Ogk$vectors
       Lambda.Ogk <- Eigen.R.Ogk$values
       A.Ogk <- D.Ogk%*% E.Ogk
       pre.Gamma.Ogk <- x%*%solve(t(A.Ogk))
       Gamma.Ogk <- diag(apply(pre.Gamma.Ogk, 2, get.Q.adj, phi.est=0)^2)
       Ogk.est <- A.Ogk %*% Gamma.Ogk %*% t(A.Ogk)
       return(Ogk.est)
    }

    # function which gives the matrix of the slopes between each pair of observations in a window sample, 
    # required for fast update of RM slope estimation
    get.B.t <- function(y){
        n <- length(y)
        j <- 1:n
        B.t <- matrix(NA,n,n)
        for (i in j) {
            B.t[j[j<i],i] <- B.t[i,j[j<i]] <- (y[i] - y[j[j<i]]) / (i - j[j<i])
        }
        return(B.t)
    }
    
    # function which gives the RM slope estimation based on the output of function "get.B.t"
    get.beta.RM <- function(B.t){
        beta.i <- apply(B.t,1,median, na.rm=T)
        beta.RM <- median(beta.i, na.rm=T)
        return(beta.RM)
    }
    # function which gives the RM level estimation based on the RM slope estimation given by the function "get.beta.RM"
    get.mu.RM <- function(y, beta.RM){
        n <- length(y)
        j <- 1:n
        mu.RM <- median(y - (beta.RM*(j-n)), na.rm=T)
        return(mu.RM)
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

    # obtain the scarm test statistic based on the outputs "B" and "noise.sd" of the functions "get.B.t"  and "get.Q.adj"
    scarm.test.statistic <- function(y, B, ell, r, noise.sd, phi.est){
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
        var.t <- noise.sd^2 * (v.left+v.right)
        test.statistic <- d.t/(sqrt(var.t))
        return(list(test.statistic=test.statistic, beta.left=beta.left, beta.right=beta.right))
    }
    
    # model for the correlation of two RM-slopes; depends on cross-correlation of the errors and the window widths
    corr.RM.slope.model <- function(ratio,rho){
        out <- (rho*0.08986669) + ((rho^3)*0.01312361) - (rho*ratio*1.17617121) + (rho*(ratio^2)*3.57618663) - (rho*(ratio^3)*1.72785517) - ((rho^3)*ratio*0.18067041)+ ((rho^3)*(ratio^2)*0.58289707) - ((rho^3)*(ratio^3)*0.19580252)
        return(out)
    }
    
    # model for the variance of the RM-slope for a linear signal with Gaussian error in a window sample of length "width"
    var.RM.slope.model <- function(width){
        out <- 4.77e-07 + (17.71*(1/width^3))
        return(out)
    }
    
    # function for block building based on the comparison statistics in "S.est"
    build.blocks <- function(S){
        S <- abs(S)
        S[is.na(S)] <- c.bound+1
        K <- 0.5 + sqrt(2*length(S)+0.25)
        j <- 1 
        R <- 1:K
        B <- list()
        combs <- combn(R,2)
        repeat{
            if(all(S >= c.bound)){
                if(length(R)>0){
                    for(i in 1:length(R)){
                        B[[j+i-1]] <- R[i]
                    }
                }
                K.t <- length(B)
                break
            }
            if(all(S < c.bound)){
                B[[j]] <- R
                K.t <- length(B)
                break
            }
            if(any(S < c.bound)){
                pos <- which.min(S)
                B[[j]] <- combs[,pos]
                R <- R[is.na(match(R,B[[j]]))]
                S[pos] <- Inf
                combs[,pos] <- 0
            }
            repeat{
                if(all(S >= c.bound)){
                    K.t <- length(B)
                    break
                }
                aligned <- integer()
                aligned.pos <- integer()
                for(i in 1:length(R)){
                    k <- R[i]
                    check1 <- which(apply(combs,2,function(x){any(!is.na(match(x,k))) & any(!is.na(match(x,B[[j]])))}))
                    check2 <- which(apply(combs,2,function(x){any(!is.na(match(x,k))) & any(!is.na(match(x,B[[j]])))}) & (S < c.bound))
                    if(length(check1)==length(check2) & length(check1)>0 & length(check2)>0 ){
                        aligned <- c(aligned,k)
                        aligned.pos <- c(aligned.pos,check1)
                    }
                }
                if(length(aligned)==0){
                    for(i in 1:dim(combs)[2]){
                        if(any(B[[j]]==combs[1,i]) | any(B[[j]]==combs[2,i])){
                            combs[,i] <- 0
                            S[i] <- Inf
                        }
                    }   
                    break
                }
                best.aligned.pos <- min(which(S==min(S[aligned.pos])))
                if(any(R==combs[1,best.aligned.pos])){
                    k <- combs[1,best.aligned.pos]
                } else {
                    k <- combs[2,best.aligned.pos]
                }
                B[[j]] <- sort(c(B[[j]],k))
                R <- R[-which(R==k)]
                for(i in 1:length(B[[j]])){
                    k <- B[[j]][i]
                    remove.pos <- apply(combs,2,function(x){any(x==k)})
                    combs[,remove.pos] <- 0
                    S[remove.pos] <- Inf
                }
            }
            if(length(R)==0){
                K.t <- length(B)
                break
            }
            j <- j+1
        }
        blocks <- rep(0,K)
        for(i in 1:K.t){
            blocks[B[[i]]] <- i
        }
        return(blocks)
    }

    # function to allocate a color to each block
    block.col <- function(x){ 
        cols <- x
        color <- 1
        for(i in 1:length(unique(x))){
            if(length(which(x==i))==1){
                cols[x==i]<- 0
            } else {
                cols[x==i] <- color
                color <- color + 1
            }
        }
        return(cols)
    }
    
    # function to sort the colors which are allocated to the blocks
    sort.colors <- function(x){
        if(all(is.na(x)) | all(x==0)){
            return(x)
        } else {
            n <- length(x)
            not.zero <- which(x!=0)
            y <- x[not.zero]
            K.t <- max(x)
            if(K.t>1){
                y <- -y
                m <- length(y)
                z <- c()
                i <- 1
                repeat{
                    block.i <- which(y==y[1])
                    z <- c(z,rep(i,length(block.i)))
                    y <- y[-block.i]
                    if(length(y)==0) break
                    i <- i+1
                }
                x[not.zero] <- z
            }
            return(x)
        }
    }
    
    #################
    # Main function #
    #################
    calc.steps <- dim.series*(dim.series-1)/2
    combs <- combn(1:dim.series,2)
    
    # empty vectors, matrices and lists for results
    scarm.width.mat <- scarm.signal.mat <- slope.mat <- noise.sd.mat <- scarm.statistic <- scarm.critval <- acf.lag.one.mat <- matrix(NA, ncol=dim.series, nrow=length.series)
    blocks <- block.colors <- matrix(NA, ncol=dim.series, nrow=length.series)
    signal.mat <- width.mat <- matrix(NA, ncol=dim.series, nrow=length.series)
    width.mat[n.min-1,] <- n.min-1
    S.est <- S.est2 <- matrix(NA, ncol=calc.steps, nrow=length.series)
    Cov.est <- Corr.est <- NULL
    K.t <- rep(NA, length.series)
    width.mat[n.min-1,] <- n.min-1
    
    # matrix of the pairwise slopes, obtained by "get.B.t"; required for the very first signal estimation
    B.matrices <- NULL
    for(k in 1:dim.series){
        B.matrices[[k]] <- get.B.t(time.series[1:n.min,k])
    }

    # begin loop over time points i
    for(i in n.min:length.series){
    
        # 1.) apply SCARM for each univariate component of the multivariate time series
        for(k in 1:dim.series){
            
            # compute autocorrelation at lag 1 for each variable at time point i
            phi.est.sample.data <- time.series[(max(i-n.max+1,1)):i,k]
            phi.est.sample.signal <- signal.mat[(max(i-n.max+1,1)):i,k]
            phi.est.sample <- na.omit(phi.est.sample.data-phi.est.sample.signal)
            if(length(phi.est.sample)>=n.min){
                phi.est <- acf(phi.est.sample, plot=F)$acf[2]
                if(is.na(phi.est)){
                    phi.est <- 0
                }
            } else {
                phi.est <- 0
            }
            if(is.na(slope.mat[i-1,k])){
                n <- min(width.mat[i-1,k]+1, ell.min+r)
            } else {
                n <- min(width.mat[i-1,k]+1,n.max)
            }
            
            y <- time.series[(i-n+1):i,k]
            B <- B.matrices[[k]]
            if(all(is.na(y[(n-n.min+1):n])) || length(which(!is.na(y[(n-(min(r,n))+1):n]))) < n.min){
                # there are not enough recent observations!
                mu.est <- NA
                beta.est <- NA
                noise.sd.est <- NA
                # Update step
                m <- min(n,ell.min+r)
                #adapted.width[i,k] <- m
                if(i != length.series){
                    b.new <- (time.series[i+1,k] - time.series[i+c((-m + 1):0),k]) / (i+1 - (i+c((-m + 1):0)))
                    B <- B[(ncol(B)-m+1):ncol(B),(ncol(B)-m+1):ncol(B)]
                    B <- cbind(rbind(B, b.new), c(b.new, NA))
                    #y     <- time.series[(i - m + 1):(i + 1)]
                    #n     <- m + 1
                }
            } else {
                if(n < ell.min+r){
                # n is not large enough for test application!
                    noise.sd.est <- get.Q.adj(y, phi.est=phi.est)
                    beta.est <- get.beta.RM(B)
                    mu.est <- get.mu.RM(y=y, beta.RM=beta.est)
                    if(rtr==TRUE){
                        rtr.sample <- time.series[(i-n.min+1):i,k]
                        min.rtr.sample <- min(rtr.sample, na.rm=T)
                        max.rtr.sample <- max(rtr.sample, na.rm=T)
                        mu.est <- max(min.rtr.sample, min(mu.est, max.rtr.sample))
                    }
                    # Update step
                    if(i != length.series){
                        if(n == n.max){ 
                        # is window width maximal?
                            b.new <- (time.series[i+1,k] - time.series[i+c((-n + 2):0),k]) / (i+1 - (i+c((-n + 2):0)))
                            B <- cbind(rbind(B[-1, -1], b.new), c(b.new,NA))
                        } else {
                            b.new <- (time.series[i+1,k] - time.series[i+c((-n + 1):0),k]) / (i+1 - (i+c((-n + 1):0)))
                            B <- cbind(rbind(B, b.new), c(b.new, NA))
                        }
                    }       
                } else {
                  # is n large enough for test application?
                    ell <- n - r
                    if(length(which(!is.na(y[1:ell])))<ell.min/2 || length(which(!is.na(y[(ell+1):n])))<r/2){
                    # there are not enough observations for testing!
                        B <- B[(ncol(B)-ell.min-r+1):ncol(B),(ncol(B)-ell.min-r+1):ncol(B)]
                        noise.sd.est <- get.Q.adj(y, phi.est=phi.est)
                        beta.est <- get.beta.RM(B)
                        mu.est <- get.mu.RM(beta.RM=beta.est, y=y[(n-ell.min-r):n])
                        n <- ell.min+r
                        if(rtr==TRUE){
                            rtr.sample <- time.series[(i-n.min+1):i,k]
                            min.rtr.sample <- min(rtr.sample, na.rm=T)
                            max.rtr.sample <- max(rtr.sample, na.rm=T)
                            mu.est <- max(min.rtr.sample, min(mu.est, max.rtr.sample))
                        }
                        if(i != length.series){
                            b.new <- (time.series[i+1,k] - time.series[i+c((-n + 1):0),k]) / (i+1 - (i+c((-n + 1):0)))
                            B <- cbind(rbind(B, b.new), c(b.new,NA))
                        }
                    } else {
                    # are there enough observations for testing? 
                        noise.sd.est <- get.Q.adj(y, phi.est=phi.est)
                        scarm.critval[i,k] <- get.critval(ell=ell, r=r, sign.level=sign.level)
                        scarm.test <- scarm.test.statistic(y=y, B=B, ell=ell, r=r, noise.sd=noise.sd.est, phi.est=phi.est)
                        scarm.statistic[i,k] <- scarm.test$test.statistic
                        slope.diff <- scarm.test$beta.left - scarm.test$beta.right
                        if(abs(scarm.statistic[i,k]) > scarm.critval[i,k]){
                        # does test reject the null hypothesis?
                            B <- B[(ncol(B)-n.min+1):ncol(B),(ncol(B)-n.min+1):ncol(B)]
                            beta.est <- get.beta.RM(B)
                            mu.est <- get.mu.RM(beta.RM=beta.est, y=y[(n-n.min+1):n])
                            n <- n.min
                            if(rtr==TRUE){
                                rtr.sample <- time.series[(i-n.min+1):i,k]
                                min.rtr.sample <- min(rtr.sample, na.rm=T)
                                max.rtr.sample <- max(rtr.sample, na.rm=T)
                                mu.est <- max(min.rtr.sample, min(mu.est, max.rtr.sample))
                            }
                            if(i != length.series) {
                                b.new <- (time.series[i+1,k] - time.series[i+c((-n + 1):0),k]) / (i+1 - (i+c((-n + 1):0)))
                                B <- cbind(rbind(B, b.new), c(b.new,NA))
                            }
                        } else {
                        # difference of RM slopes is NOT too large => do not decrease time window
                            beta.est <- get.beta.RM(B)
                            mu.est <- get.mu.RM(beta.RM=beta.est, y=y)
                            if(i != length.series){
                                if(n == n.max){ 
                                # is window width maximal?
                                    b.new <- (time.series[i+1,k] - time.series[i+c((-n + 2):0),k]) / (i+1 - (i+c((-n + 2):0)))
                                    B <- cbind(rbind(B[-1, -1], b.new), c(b.new,NA))
                                } else {
                                    b.new <- (time.series[i+1,k] - time.series[i+c((-n + 1):0),k]) / (i+1 - (i+c((-n + 1):0)))
                                    B <- cbind(rbind(B, b.new), c(b.new, NA))
                                }
                            }
                            if(rtr==TRUE){
                                rtr.sample <- time.series[(i-n.min+1):i,k]
                                min.rtr.sample <- min(rtr.sample, na.rm=T)
                                max.rtr.sample <- max(rtr.sample, na.rm=T)
                                mu.est <- max(min.rtr.sample, min(mu.est, max.rtr.sample))
                            }
                        }
                    }
                }
            }
            # outputs of the SCARM
            B.matrices[[k]] <- B
            scarm.width.mat[i,k] <- n
            scarm.signal.mat[i,k] <- mu.est
            slope.mat[i,k] <- beta.est
            noise.sd.mat[i,k] <- noise.sd.est
            acf.lag.one.mat[i,k] <- phi.est
            width.mat[i,] <- scarm.width.mat[i,]
        }
              
        # treatment of missing values
        non.NA.pos <- which(!is.na(slope.mat[i,]))
        dim.i <- length(non.NA.pos)
        if(dim.i==1){
            signal.mat[i,non.NA.pos] <- scarm.signal.mat[i,non.NA.pos]
            width.mat[i,non.NA.pos] <- scarm.width.mat[i,non.NA.pos]
        }
        if(dim.i > 1){
            replace.pos <- apply(combs, 2 , function(x){any(x[1]==non.NA.pos) & any(x[2]==non.NA.pos)})
            scarm.width.i <- scarm.width.mat[i,non.NA.pos]
            slope.i <- slope.mat[i,non.NA.pos]
            noise.sd.i <- noise.sd.mat[i,non.NA.pos]
            combs.i <- combn(1:dim.i,2)
            calc.steps.i <- dim.i*(dim.i-1)/2
            korr.beta <- corr.ogk.est <- rep(NA,calc.steps.i)
            
            # 2.) obtain the comparison statistics in "S.est" for each time point i; 
            # use SCARM-outputs "slope.mat", "noise.sd.mat" and "scarm.width.mat" from step 1.)
            
            # compute the ration "width.ratio" for each pair of components (smaller width is divided by larger width)
            max.w <- apply(combs.i, 2, function(x){apply(cbind(scarm.width.i[x[1]],scarm.width.i[x[2]]),1,max)})
            min.w <- apply(combs.i, 2, function(x){apply(cbind(scarm.width.i[x[1]],scarm.width.i[x[2]]),1,min)})
            width.ratio <- min.w/max.w
            # estimate the matrix "Cov.est" and "Corr.est" of the error covariances and error correlations
            Cov.est[[i]] <- Corr.est[[i]] <- matrix(NA,dim.series,dim.series)
            n <- max(scarm.width.i)
            sample.matrix <- time.series[(i-n+1):i,non.NA.pos]
            if(any(is.na(sample.matrix))){
                noise.sd.est <- as.vector(noise.sd.mat[i,non.NA.pos])
                RM.slopes <- as.vector(slope.mat[i,non.NA.pos])
                RM.levels <- as.vector(scarm.signal.mat[i,non.NA.pos])
                RM.lines <- matrix(NA, nrow=n, ncol=length(non.NA.pos))
                for(k in 1:length(non.NA.pos)){
                    RM.lines[,k] <- (RM.levels[k] + (RM.slopes[k]*((1-n):0)))
                    na.pos <- which(is.na(sample.matrix[,k]))
                    #set.seed(i*k)
                    sample.matrix[na.pos,k] <- RM.lines[na.pos,k]# + rnorm(n=length(na.pos), sd=noise.sd.est[k])
                }  
            }
            Ogk.est <- get.Ogk(sample.matrix, noise.sd.est.s=noise.sd.mat[i,non.NA.pos])
            Cov.est[[i]][non.NA.pos,non.NA.pos] <- Ogk.est
            Corr.est[[i]][non.NA.pos,non.NA.pos] <- cov2cor(Ogk.est)
            # estimate the correlation "korr.beta" between each pair of RM-slopes
            # by using the estimated error correlations and the model "corr.RM.slope.model"
            for(m in 1:calc.steps.i){  
                korr.beta[m] <- corr.RM.slope.model(ratio=width.ratio[m], rho=cov2cor(Ogk.est)[combs.i[1,m],combs.i[2,m]])
            }
            # estimate the variances of the RM-slopes ("var.slopes")
            # by using the estimated error variance and the model "var.RM.slope.model"
            var.slopes <- apply(matrix(1:dim.i), 1, function(x){noise.sd.i[x]^2 * var.RM.slope.model(scarm.width.i[x])})
            # estimate the variance of the difference between each pair of RM-slopes
            var.slope.diffs <- apply(combs.i, 2, function(x){var.slopes[x[1]] + var.slopes[x[2]] - 2*(korr.beta[which(combs.i[1,]==x[1]&combs.i[2,]==x[2])]*sqrt(var.slopes[x[1]])*sqrt(var.slopes[x[2]]))})
            # estimate the comparison statistics "S.est"
            slope.diffs <- apply(combs.i, 2, function(x){slope.i[x[1]] - slope.i[x[2]]})
            S.est.i <- slope.diffs/sqrt(var.slope.diffs)
            S.est[i,replace.pos] <- S.est.i
            if(r.bound > 0){
                r.bound.Index.no <- width.ratio < r.bound
                r.bound.Index.yes <- width.ratio >= r.bound
                S.est.i2 <- (c.bound+1)*r.bound.Index.no + r.bound.Index.yes*S.est.i            
            } else {
                S.est.i2 <- S.est.i
            }
            S.est2[i,replace.pos] <- S.est.i2
                    
            # 3.) build blocks based on "S.est"
            blocks[i,] <- build.blocks(S.est2[i,])
            block.colors[i,] <- block.col(blocks[i,])
    
            # 4.) if a block contains two or more components, the multivariate TRM-LS is applied for signal estimation
            # using the estimated error covariance matrix "Cov.est" from step 2.)
        
            # blocks and number of blocks at time point i
            B <- blocks[i,]
            K.t[i] <- length(unique(B))
            # loop over all blocks
            for(j in 1:K.t[i]){
                block.comps <- which(B==j)
                block.size <- length(block.comps)
                if(block.size==1){
                    signal.mat[i,block.comps] <- scarm.signal.mat[i,block.comps]
                    width.mat[i,block.comps] <- scarm.width.mat[i,block.comps]
                }
                if(block.size>1){
                    # obtain the RM-residuals within the window of width "block.width"
                    indiv.widths <- as.vector(scarm.width.mat[i,block.comps]) 
                    block.width <- min(indiv.widths)
                    RM.slopes <- as.vector(slope.mat[i,block.comps])
                    RM.levels <- as.vector(scarm.signal.mat[i,block.comps])
                    RM.lines <- matrix(NA, nrow=block.width, ncol=block.size)
                    for(k in 1:block.size){
                        n <- indiv.widths[k]
                        RM.lines[,k] <- (RM.levels[k] + (RM.slopes[k]*((1-n):0)))[(n-block.width+1):n]
                    }
                    block.sample <- time.series[(i-block.width+1):i,block.comps]
                    block.sample[is.na(block.sample)] <- RM.lines[is.na(block.sample)]
                    res <- block.sample - RM.lines
                    # detect outlying observation vectors
                    Cov.est.j <- Cov.est[[i]][block.comps,block.comps] 
                    mahala.dist <- mahalanobis(res, center=rep(0,block.size), cov=Cov.est.j, tol=2.220446e-16)
                    cutoff <- qchisq(p=0.99, df=block.size)*median(mahala.dist)/qchisq(p=0.5, df=block.size)
                    I.j <- which(mahala.dist <= cutoff)
                    if(length(I.j) < block.width/2){
                        I.j <- 1:block.width      
                    }
                    # LS-Regression based on the trimmed window sample
                    design.points <- 1:block.width
                    LS.cov    <- cov(cbind(design.points,block.sample)[I.j,])
                    LS.slopes   <- LS.cov[1,2:(block.size+1)] / LS.cov[1,1]
                    LS.levels  <- apply(block.sample[I.j,], 2, mean) - (LS.slopes * mean(design.points))
                    LS.signals <- LS.levels + (LS.slopes*block.width)
                    # apply restrict-to-range-rule
                    if(rtr==TRUE){
                        for(k in 1:block.size){
                            LS.signal <- LS.signals[k]
                            rtr.sample <- time.series[(i-n.min+1):i,block.comps[k]] 
                            min.rtr <- min(rtr.sample, na.rm=T)
                            max.rtr <- max(rtr.sample, na.rm=T)
                            LS.signals[k] <- min(max(LS.signal, min.rtr),max.rtr)
                        }
                    }
                    # store TRM-LS signal estimations and window widths
                    signal.mat[i,block.comps] <- LS.signals
                    width.mat[i,block.comps] <- block.width
                }
            }
            # the block building may induce that some individual window widths are decreased;
            # in this case the matrices of the pairwise slopes have to be adjusted
            if(any(width.mat[i,]!=scarm.width.mat[i,])){
                new.widths <- which(width.mat[i,]!=scarm.width.mat[i,])
                for(k in 1:length(new.widths)){
                    B <- B.matrices[[new.widths[k]]]
                    old.dim.B <- dim(B)[1]
                    new.dim.B <- width.mat[i,new.widths[k]]+1
                    B.matrices[[new.widths[k]]] <- B[(old.dim.B-new.dim.B+1):old.dim.B,(old.dim.B-new.dim.B+1):old.dim.B]
                }       
            }   
        }
    }   # end of loop over time points i
    width.mat[n.min-1,] <- NA
    block.colors <- t(apply(block.colors, 1, sort.colors))

    ###########
    # Outputs #
    ###########
    result <- list(signal.est=signal.mat,
                   slope.est=slope.mat,
                   adapted.width=width.mat,
                   noise.sd.est=noise.sd.mat,
                   #
                   scarm.signal.est=scarm.signal.mat,
                   scarm.width=scarm.width.mat,
                   scarm.statistic=scarm.statistic,
                   scarm.critval=scarm.critval,
                   #
                   ssm.statistic=S.est,
                   blocks=blocks,
                   acf.lag.one=acf.lag.one.mat,
                   #
                   dim.series=dim.series,
                   length.series=length.series,
                   S.est2=S.est2,             
                   block.colors=block.colors,      
                   # all (maybe adjusted) inputs are given as outputs
                   time.series=time.series,
                   right.width=r,
                   min.left.width=ell.min,
                   min.width=n.min,
                   max.width=n.max,
                   sign.level=sign.level,
                   bound.noise.sd=bound.noise.sd,
                   rtr=rtr,
                   autocorrelations=autocorrelations,                   
                   c.bound=c.bound,
                   r.bound=r.bound
                   )             
    return(structure(result, class = "mscarm.filter"))
}
   
   
##################
# Default output #
##################
print.mscarm.filter <- function(x, ...){
    N <- dim(x$signal.est)[1]
    cat('$signal.est \n')
    if(N <= 100){
        print(x$signal.est, ...)
    } else {
        print(x$signal.est[1:100,])
        cat('Only the first 100 signal estimations are printed.\n')
    }
}
   
########
# Plot #
########

plot.mscarm.filter <- function(x, info=FALSE,...){
    if(info==TRUE){
        par(mfrow=c(2,1))
    } else {
        par(mfrow=c(1,1))
    }
    y.lim.min <- min(x$time.series,na.rm=T)
    y.lim.max <- max(x$time.series,na.rm=T)
    y.lim.diff <- y.lim.max-y.lim.min
    plot(x$time.series[,1], ylim=c(y.lim.min-(y.lim.diff/5),y.lim.max+(y.lim.diff/5)), ylab="", main="Multivariate time series and mSCARM signal extraction", col="grey", xlab="t", type="l")
    points(x$time.series[,1], pch=20, cex=0.5, col="grey")
    lines(x$signal.est[,1])
    legend("topleft", c("Data", "Signal extraction"), col=c("grey","black"), lty=c(1,1), lwd=c(1,1), bty="n")
    #if(info==TRUE){
    #    legend("topright","Detected signal changes", col=1, pch=4, lty=NA, lwd=2, bty="n")
    #}
    for(i in 1:x$dim.series){
        text(paste("V",i,sep=""), x=1, y=na.omit(x$time.series[,i])[1], pos=2)
    }
    for(j in 2:x$dim.series){
        lines(x$time.series[,j], col="grey")
        points(x$time.series[,j], pch=20, cex=0.5, col="grey")
        lines(x$signal.est[,j])
    }
    if(info==TRUE){
        #for(i in 1:x$dim.series){
        #    breaks.i <- which(abs(x$scarm.statistic)[,i] > x$scarm.critval[,i]) - floor(x$min.width/3)
        #    points(x=breaks.i, y=x$signal.est[breaks.i,i], pch=4, cex=1, lwd=2)
        #}
        comparisons <- (x$dim.series)*(x$dim.series-1)/2
        col.names <- c()
        j <- 1
        for(i in seq(1,comparisons*2,2)){
            col.names[j] <- paste("V",combn(x$dim.series,2)[i],"-V",combn(x$dim.series,2)[i+1], sep="")
            j <- j+1
        }
        plot(rep(NA,x$length.series), type="l", ylim=c(1,comparisons+1), xlab="t", ylab="", yaxt="n", 
        main=paste("Pairwise slope comparison between the variables V1,...,V",x$dim.series, sep=""))
        axis(2, at=1:comparisons, labels=col.names, las=1)
        for(i in 1:comparisons){
            lines(x=1:(x$length.series),y=rep(i,x$length.series), col="lightgrey")
            crosses <- rep(NA,x$length.series)
            crosses[which(abs(x$S.est2[,i]) <= x$c.bound)] <- i
            points(crosses, pch=3)
        }
        legend("topleft", "Similar slope detected", lty=NA, lwd=NA, pch=3, bty="n", col=1)
    }
}
                
                   
                   
