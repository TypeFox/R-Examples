Bnk.func = function(x0, x, h, kernel.type="4th"){
    # the kernel weight function Bnk(x0, x), where x0 is a scalar, and x is a vector
    # returns a vector
    # h is the bandwidth
    if(!is.vector(x0)) x0 <- as.vector(x0)

    xx <- (x-x0)/h 
    if(kernel.type == "4th"){ 
        xx[abs(xx) >= 1] <- 1
        w <- 15 * (1 - xx^{2})^{2/16}  #biquadratic kernel 
    }
    w <- w/sum(w)
    return(w)
}

Kernel.func <- function(U0, U, h, kernel.type="4th"){
    # U: n*k matrix
    # U0: 1-k matrix
    # return: K((U-U0)/h)
    if(!is.vector(U0)) U0 <- as.vector(U0)
    n = nrow(U)
    if(kernel.type=="4th"){
        tt = rbind(U, U0)
        tmp = apply(tt, 2, function(x) {
                Bnk.func(x0 = x[n+1], x = x[1:n], h = h, kernel.type = kernel.type)
            })
        tmp = apply(tmp, 1, prod)
        tmp = tmp/sum(tmp)
    }
    return(tmp)    
}

tauhat.func <- function(y0, x0, y, x, delta,h, kernel.type = "4th"){
# tau0(y0, x0) = F(T<y0|x0); so y0 is the C_i, and x0 is the xi in the paper
# x0: k-dimensional covariate vector
# y: n-vector of observed survival time = T^C
# x: k-dimensional covariate matrix
# delta: the censoring indicator function
# h: bandwidth parameter
    n <- length(y)
    w <- rep(0,n)
    ###kernel weights#########################
    p <- qr(x)$rank

    if(p >1) Bn = Kernel.func(x0, x, h, kernel.type)
    else if (p == 1) Bn = Bnk.func(x0, x, h, kernel.type)

    if (y0<max(y)) {
        # sort the data y, and the delta, Bn correspondingly to the order of sorted y
        y2 = sort(y)
        Order = order(y) # so z[Order] = z2
        Bn2 = Bn[Order]
        delta2 = delta[Order]
        eta = which(delta2==1 & y2<=y0) # the index of those observations satisfying delta2==1 & z2<=y0
        Bn3 = Bn2[n:1]  # change the order of Bn2, make the first obs of Bn2 to be the last of Bn3
        tmp = 1 - Bn2 /cumsum(Bn3)[n:1]  
        out = 1 - prod(tmp[eta], na.rm=T) # na.rm=T, as some of those tmp=NA as the denom =0
    } 
    else out<-1 
    return(out)
}

LCRQ <- function(y, x, delta, tau, h, kernel.type = "4th"){
# Locally weighted censored quantile regression method
# x is a design matrix
# y is the observed survival time = min(T, C)
# delta is the censoring indicator function with 1 standing for uncensored, and 0 censored
# tau is the quantile level of interest
# h is the handwidth used in calculating tauhat for each censored observation
    if(!is.matrix(x)) x <- as.matrix(x)

    n <- length(y)
    ind <- which(delta == 0)
    w <- rep(1, n) # the weight vector

    if(length(ind) >= 1){
        for(i in 1:length(ind)){
            x0 = x[ind[i], , drop = FALSE]
            y0 = y[ind[i]]
            tau.star <- tauhat.func(y0, x0, y, x, delta, h =h, kernel.type = "4th")
            if (tau > tau.star) w[ind[i]] <- (tau - tau.star) / (1-tau.star)
            else w[ind[i]] <- 0
        }
        # pseudo observations
        ind2 <- which(w != 1)
        y.pse <- rep(max(y)+100, length(ind2))
        x.pse <- x[ind2, , drop = FALSE]

        yy <- c(y, y.pse)
        xx <- rbind(x, x.pse)
        ww <- c(w, 1-w[ind2])
    }
    else{
        yy = y
        xx = x
        ww = w
    }
    
    rq1 = rq(yy~xx, weights=ww, tau=tau)
    result<-rq1$coeff
    result    
}
