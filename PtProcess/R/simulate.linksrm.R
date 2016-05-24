simulate.linksrm <- function(object, nsim=1, seed=NULL, max.rate=NA,
                             stop.condition=NULL, ...){
    data <- object$data
    params <- object$params
    gparams <- eval(object$gmap)
    mparams <- eval(object$mmap)
    gif <- object$gif
    TT <- object$TT
    if (!is.null(attr(gif, "regions"))) n <- eval(attr(gif, "regions"))
    else stop("gif needs a regions attribute")
    if (!is.null(data)) {
        use <- (data[, "time"] < TT[1])
        if (sum(use) == 0) 
            data <- NULL
        else data <- data[use, c("time", "magnitude", "region")]
    }
    set.seed(seed)
    ti <- TT[1]
    repeat {
        Rmax <- 0
        rate <- sum(gif(data = data, evalpts = cbind(time = rep(ti, 
            n), region = 1:n), params = gparams, tplus = TRUE))
        tmax <- ti
        while (rate > Rmax) {
            ti <- tmax
            tmax <- ti + qexp(0.4, rate = rate)
            Rmax <- sum(gif(data = data, evalpts = cbind(time = rep(tmax, 
                n), region = 1:n), params = gparams))
            tau <- rexp(1, rate = Rmax)
            rate <- sum(gif(data = data, evalpts = cbind(time = rep(ti + 
                tau, n), region = 1:n), params = gparams))
        }
        ti <- ti + tau
        if (ti > TT[2]) 
            break
        if (runif(1, 0, 1) <= rate/Rmax) {
            #    accept simulated point at ti
            #    now generate accompanying marks
            newevent <- list()
            newevent$time <- ti
            #    simulate the region
            bound <- gif(data, evalpts = cbind(time = rep(ti, 
                n), region = 1:n), gparams)
            bound <- cumsum(bound)/sum(bound)
            newevent$region <- sum(bound < runif(1, 0, 1)) + 1
            #    simulate other marks
            newevent <- c(newevent, object$marks[[2]](ti, data, mparams))
            newevent <- as.data.frame(newevent)
            data <- rbind(data, newevent)
            if (!is.null(stop.condition))
                if (stop.condition(data)) break
        }
    }
    object$data <- data
    return(object)
}
