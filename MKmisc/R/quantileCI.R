## Confidence Intervals for quantiles
quantileCI <- function(x, prob = 0.5, conf.level = 0.95, method = "exact", 
                       minLength = FALSE, na.rm = FALSE){
    cl <- match.call()
    if(!is.na(pmatch(method, "exact")))
        method <- "exact"

    METHODS <- c("exact", "asymptotic")
    method <- pmatch(method, METHODS)

    if(is.na(method))
        stop("invalid method")

    if(method == -1)
        stop("ambiguous method")

    if(length(x) <= 1)
        stop("'x' has to be of at least length 2")
    if(length(conf.level) != 1)
        stop("'conf.level' has to be of length 1 (confidence level)")
    if(conf.level < 0.5 | conf.level > 1)
        stop("'conf.level' has to be in [0.5, 1]")

    alpha <- 1 - conf.level
    z <- qnorm(1-alpha/2)
    if(na.rm) x <- x[!is.na(x)]
    n <- length(x)
    est <- median(x)
    xs <- sort(x)
    
    if(method == 1){ # exact
        CI.mat <- matrix(NA, ncol = 2, nrow = n-1)
        pcov.vec <- numeric(n-1)
        for(i in 1:(n-1)){
          for(j in (i+1):n){
            pcov <- pbinom(j-1, size = n, prob = prob)-pbinom(i-1, size = n, prob = prob)
            if(pcov > conf.level){
              pcov.vec[i] <- pcov
              CI.mat[i,] <- c(xs[i], xs[j])
              break
            }
          }
        }
        if(all(pcov.vec == 0)){
          CI <- c(xs[1], xs[n])
          attr(CI, "confidence level") <- 1
        }else{
          CI.mat <- CI.mat[pcov.vec > 0,,drop = FALSE]
          pcov.vec <- pcov.vec[pcov.vec > 0]
          pcov.min <- min(pcov.vec)
          CI <- CI.mat[pcov.vec == pcov.min,]
          if(minLength){
            CI <- CI[which.min(diff(t(CI))),]
          }
          attr(CI, "(exact) confidence level") <- pcov.min
        }
    }
    if(method == 2){ # approx
        prob.sd <- sqrt(n*prob*(1-prob))
        k.lo <- max(1, floor(n*prob - z*prob.sd))
        k.up <- min(n, ceiling(n*prob + z*prob.sd))
        CI <- c(xs[k.lo], xs[k.up])
        attr(CI, "(asymptotic) confidence level") <- conf.level
    }

    list("call" = cl, "estimate" = est, "CI" = CI)
}

medianCI <- function(x, conf.level = 0.95, method = "exact", minLength = FALSE, na.rm = FALSE){
    cl <- match.call()
    res <- quantileCI(x, prob = 0.5, conf.level = conf.level, method = method, 
                      minLength = minLength, na.rm = na.rm)
    res$call <- cl
    res
}
