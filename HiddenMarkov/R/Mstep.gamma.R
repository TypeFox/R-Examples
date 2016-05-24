"Mstep.gamma" <-
function (x, cond, pm, pn, maxiter = 200) 
{
    m <- ncol(cond$u)
    nms <- sort(names(pm))
    if (all(nms == c("rate", "shape"))) {
        denom <- apply(cond$u, MARGIN = 2, FUN = sum)
        y1 <- as.numeric(matrix(x, nrow = 1) %*% cond$u)/denom
        y2 <- as.numeric(matrix(log(x), nrow = 1) %*% cond$u)/denom
        rate <- pm$rate
        shape <- pm$shape
        dLL <- matrix(NA, ncol = 1, nrow = 2)
        info <- matrix(NA, ncol = 2, nrow = 2)
        for (j in 1:m) {
            p <- matrix(c(rate[j], shape[j]), ncol = 1)
            for (iter in 1:maxiter) {
                dLL[1, 1] <- p[2, 1]/p[1, 1] - y1[j]
                dLL[2, 1] <- y2[j] - digamma(p[2, 1]) + log(p[1, 
                  1])
                info[1, 1] <- -p[2, 1]/p[1, 1]^2
                info[1, 2] <- 1/p[1, 1]
                info[2, 1] <- 1/p[1, 1]
                info[2, 2] <- -trigamma(p[2, 1])
                info <- solve(info)
                incr <- info %*% dLL
                p <- p - incr
                if (any(p <= 0))
                    stop("Parameter estimates out of bounds")
                if (all(abs(incr) < 1e-05)) 
                    break
            }
            if (iter == maxiter)
                warning("Maximum iterations reached")
            rate[j] <- p[1, 1]
            shape[j] <- p[2, 1]
        }
        return(list(rate = rate, shape = shape))
    }
    if (all(nms == "rate")) {
        rate <- as.numeric(matrix(pn$shape, nrow = 1) %*% cond$u)/as.numeric(matrix(x, 
            nrow = 1) %*% cond$u)
        return(list(rate = rate))
    }
    if (all(nms == "shape")) {
        denom <- apply(cond$u, MARGIN = 2, FUN = sum)
        y2 <- as.numeric(matrix(log(x), nrow = 1) %*% cond$u)/denom
        lograte <- as.numeric(matrix(log(pn$rate), nrow = 1) %*% 
            cond$u)/denom
        shape <- pm$shape
        for (j in 1:m) {
            for (iter in 1:maxiter) {
                dLL <- y2[j] - digamma(shape[j]) + lograte[j]
                incr <- -dLL/trigamma(shape[j])
                shape[j] <- shape[j] - incr
                if (shape[j] <= 0)
                    stop("Parameter estimates out of bounds")
                if (all(abs(incr) < 1e-05)) 
                    break
            }
            if (iter == maxiter)
                warning("Maximum iterations reached")
        }
        return(list(shape = shape))
    }
    stop("Invalid specification of parameters")
}
