CV.S=function (y, S, W = NULL, trim = 0, draw = FALSE, metric = metric.lp, ...) {
    n = ncol(S)
    isfdata <- is.fdata(y)
    if (isfdata) {
        nn<-nrow(y)
        if (is.null(W)) W<-diag(nn)
        y2 = t(y$data)
        y.est = t(S %*% y2)
        y.est <- fdata(y.est, y$argvals, y$rangeval, y$names)
        e <- (y - y.est)/(1-diag(S))
#        e$data<-sqrt(W)%*%(e$data)     
        ee <- drop(norm.fdata(e, metric = metric, ...)[, 1]^2)
        if (trim > 0) {
            e.trunc = quantile(ee, probs = (1 - trim), na.rm = TRUE, type = 4)
            ind <- ee <= e.trunc
            if (draw)   plot(y, col = (2 - ind))
            res = mean(ee[ind])
        }
        else res = mean(ee)
    }
    else {
       if (is.null(W)) W<-diag(n)
        y2 <- y
        y.est = S %*% y2
        I = diag(n)/(1 - diag(S))^2
        W = W * I
        e <- y2 - y.est
        if (trim > 0) {
            ee = t(e)
            e.trunc = quantile(abs(ee), probs = (1 - trim), na.rm = TRUE, 
                type = 4)
            l <- which(abs(ee) <= e.trunc)
            res = mean(diag(W)[l] * e[l]^2)
        }
        res = mean(diag(W) * e^2)
    }
    if (is.nan(res))    res = Inf
    return(res)
}
