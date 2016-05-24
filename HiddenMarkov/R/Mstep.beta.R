Mstep.beta <- function(x, cond, pm, pn, maxiter = 200){
    m <- ncol(cond$u)
    nms <- sort(names(pm))
    if (all(nms==c("shape1", "shape2"))){
        denom <- apply(cond$u, MARGIN=2, FUN=sum)
        y1 <- as.numeric(matrix(log(x), nrow=1) %*% cond$u)/denom
        y2 <- as.numeric(matrix(log(1-x), nrow=1) %*% cond$u)/denom
        #   use current values as initial values
        shape1 <- pm$shape1
        shape2 <- pm$shape2
        for (j in 1:m){
            p <- matrix(c(shape1[j], shape2[j]), ncol=1)
            y <- matrix(c(y1[j], y2[j]), ncol=1)
            for (iter in 1:maxiter){
                dLL <- digamma(sum(p)) - digamma(p) + y
                info <- solve(diag(-c(trigamma(p))) + trigamma(sum(p)))
                incr <- info %*% dLL
                p <- p - incr
                if (any(p <= 0))
                    stop("Parameter estimates out of bounds")
                if (all(abs(incr) < 0.00001)) break
            }
            if (iter == maxiter)
                warning("Maximum iterations reached")
            shape1[j] <- p[1,1]
            shape2[j] <- p[2,1]
        }
        return(list(shape1=shape1, shape2=shape2))
    }
    if (all(nms=="shape1")){
        denom <- apply(cond$u, MARGIN=2, FUN=sum)
        y <- as.numeric(matrix(log(x), nrow=1) %*% cond$u)/denom
        shape1 <- pm$shape1
        shape2 <- pn$shape2
        for (j in 1:m){
            uj <- as.numeric(cond$u[,j])
            dj <- as.numeric(denom[j])
            for (iter in 1:maxiter){
                D1 <- sum(digamma(shape1[j] + shape2) * uj)/dj
                D2 <- sum(trigamma(shape1[j] + shape2) * uj)/dj
                dLL <- y[j] + D1 - digamma(shape1[j])
                incr <- dLL/(D2 - trigamma(shape1[j]))
                shape1[j] <- shape1[j] - incr
                if (shape1[j] <= 0)
                    stop("Parameter estimates out of bounds")
                if (abs(incr) < 0.00001) break
            }
            if (iter == maxiter)
                warning("Maximum iterations reached")
        }
        return(list(shape1=shape1))
    }
    if (all(nms=="shape2")){
        denom <- apply(cond$u, MARGIN=2, FUN=sum)
        y <- as.numeric(matrix(log(1-x), nrow=1) %*% cond$u)/denom
        shape1 <- pn$shape1
        shape2 <- pm$shape2
        for (j in 1:m){
            uj <- as.numeric(cond$u[,j])
            dj <- as.numeric(denom[j])
            for (iter in 1:maxiter){
                D1 <- sum(digamma(shape1 + shape2[j]) * uj)/dj
                D2 <- sum(trigamma(shape1 + shape2[j]) * uj)/dj
                dLL <- y[j] + D1 - digamma(shape2[j])
                incr <- dLL/(D2 - trigamma(shape2[j]))
                shape2[j] <- shape2[j] - incr
                if (shape2[j] <= 0)
                    stop("Parameter estimates out of bounds")
                if (abs(incr) < 0.00001) break
            }
            if (iter == maxiter)
                warning("Maximum iterations reached")
        }
        return(list(shape2=shape2))
    }
    stop("Invalid specification of parameters")
}

