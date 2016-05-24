distStepSize <- function(densFn, dist,
                         param = NULL, side = c("right","left"), ...){
    set.seed(123)
    rfun <- match.fun(paste("r", densFn, sep = ""))
    side <- match.arg(side)
    n <- 50
    if(is.null(param)){
        sample <-  rfun(n, ...)
    } else {
        sample <- rfun(n, param = param)
    }
    mid <- median(sample)
    if (densFn == "skewhyp"){
        l <- list(...)
        delta <- ifelse(is.null(param), l$delta, param[2])
        nu <- ifelse(is.null(param), l$nu, param[4])
        beta <- ifelse(is.null(param), l$beta, param[3])
        if (beta > 0){
        step <- ifelse(side == "left", delta,
                       delta*abs(beta)*(nu*dist)^(-2/nu))
    }
    if (beta < 0){
        step <- ifelse(side == "right", delta,
                       delta*abs(beta)*(nu*dist)^(-2/nu))
    }
    if (isTRUE(all.equal(beta, 0))){
        step <- exp(dist/nu)
    }
        step <- c(step, mid)
    } else{
        quans <- as.vector(quantile(sample, probs = c(0.25, 0.75)))
        step <- ifelse(side == "left", mid - quans[1],
                       quans[2] - mid)
        step <- c(step, mid)
    }
    return(step)
}

