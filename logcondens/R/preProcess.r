preProcess <- function(x, xgrid = NULL){

    n <- length(x)
    if (is.null(xgrid)){
    
        # Default setting:
        xx <- sort(x)
        tmp <- c(xx[1:(n - 1)] < xx[2:n], TRUE)
        ww <- c(0, (1:n)[tmp])
        ww <- (ww[2:length(ww)] - ww[1:(length(ww) - 1)]) / n
        xx <- xx[tmp]
    }
    
    if (!is.null(xgrid)){

        if (length(xgrid) == 1){
        
        # xgrid is the distance between consecutive grid points:
        
            xx <- (floor(min(x) / xgrid):ceiling(max(x) / xgrid)) * xgrid
            } else {
                xx <- xgrid
                if (min(x) < xx[1]){xx <- c(min(x), xx)}
                if (max(x) > xx[length(xx)]){xx <- c(xx, max(x))}
            }

        m <- length(xx)
        ww <- rep(0, m)
        for (i in 1:n){
            tmp <- max(min(sum(xx <= x[i]), m - 1), 1)
            
            # tmp is an index in {1, 2, ..., m-1}
            # such that x[i] is in [xx[tmp], xx[tmp+1]].
            # Here min(..., m-1) is necessary for the
            # case of x[i] == x[mm], while max(..., 1) 
            # is for the case of x[i] < xx[1] due to
            # numerical errors.
            
            ww[tmp] <- ww[tmp] + (xx[tmp + 1] - x[i])/(xx[tmp + 1] - xx[tmp])
            ww[tmp+1] <- ww[tmp + 1] + (x[i] - xx[tmp]) / (xx[tmp + 1] - xx[tmp])
        } 
        
        ww <- ww / n
        
        # Finally, reduce xx and ww such that
        # all weights ww[i] are clearly positive:
        
        tmp <- (ww > 1e-7 / n)
        xx <- xx[tmp]
        ww <- ww[tmp]
        ww <- ww / sum(ww)
    }
    
    est.m <- sum(ww * xx)
    est.sd <- sum(ww * (xx - est.m) ^ 2)
    est.sd <- sqrt(est.sd * n / (n - 1))
    return(list("x" = xx, "w" = ww, "sig" = est.sd, "n" = n))
}













