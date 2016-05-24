"Mstep.logis" <-
function (x, cond, pm, pn, maxiter = 200) 
{
    m <- ncol(cond$u)
    nms <- sort(names(pm))
    if (all(nms == c("location", "scale"))) {
        location <- pm$location
        scale <- pm$scale
        dLL <- matrix(NA, ncol = 1, nrow = 2)
        info <- matrix(NA, ncol = 2, nrow = 2)
        for (j in 1:m) {
            p <- matrix(c(location[j], scale[j]), ncol = 1)
            for (iter in 1:maxiter) {
                denom <- sum(cond$u[,j])
                y <- 1 + exp(-(x-p[1,1])/p[2,1])
                w01 <- sum(cond$u[,j]/y)/denom
                w02 <- sum(cond$u[,j]/y^2)/denom
                w10 <- sum((x-p[1,1])*cond$u[,j])/denom
                w11 <- sum((x-p[1,1])*cond$u[,j]/y)/denom
                w12 <- sum((x-p[1,1])*cond$u[,j]/y^2)/denom
                w21 <- sum((x-p[1,1])^2*cond$u[,j]/y)/denom
                w22 <- sum((x-p[1,1])^2*cond$u[,j]/y^2)/denom
                dLL[1,1] <- (2*w01 - 1)/p[2,1]
                dLL[2,1] <- (2*w11 - w10 - p[2,1])/p[2,1]^2
                info[1,1] <- 2*(w02-w01)/p[2,1]^2
                info[1, 2] <- (1 - 2*w01)/p[2,1]^2 + 2*(w12 - w11)/p[2,1]^3
                info[2,1] <- info[1, 2]
                info[2,2] <- 1/p[2,1]^2 + 2*(w10 - 2*w11)/p[2,1]^3 + 2*(w22 - w21)/p[2,1]^4
                info <- solve(info)
                incr <- info %*% dLL
                p <- p - incr
                if (p[2,1] <= 0)
                    stop("Parameter estimates (scale) out of bounds")
                if (all(abs(incr) < 1e-05)) 
                    break
            }
            if (iter == maxiter)
                warning("Maximum iterations reached")
            location[j] <- p[1,1]
            scale[j] <- p[2,1]
        }
        return(list(location = location, scale = scale))
    }
    if (all(nms == "location")) {
        location <- pm$location
        scale <- pn$scale
        for (j in 1:m) {
            for (iter in 1:maxiter) {
                denom <- sum(cond$u[,j])
                y <- 1 + exp(-(x-location[j])/scale)
                dLL <- sum(cond$u[,j]/scale*(2/y - 1))/denom
                info <- sum(2*cond$u[,j]*(1-y)/(y*scale)^2)/denom
                incr <- dLL/info
                location[j] <- location[j] - incr
                if (all(abs(incr) < 1e-05)) 
                    break
            }
            if (iter == maxiter)
                warning("Maximum iterations reached")
        }
        return(list(location = location))
    }
    if (all(nms == "scale")) {
        scale <- pm$scale
        location <- pn$location
        for (j in 1:m) {
            for (iter in 1:maxiter) {
                denom <- sum(cond$u[,j])
                p <- scale[j]
                y <- 1 + exp(-(x-location)/p)
                w10 <- sum((x-location)*cond$u[,j])/denom
                w11 <- sum((x-location)*cond$u[,j]/y)/denom
                w21 <- sum((x-location)^2*cond$u[,j]/y)/denom
                w22 <- sum((x-location)^2*cond$u[,j]/y^2)/denom
                dLL <- (2*w11 - w10 - p)/p^2
                info <- 1/p^2 + 2*(w10 - 2*w11)/p^3 + 2*(w22 - w21)/p^4
                incr <- dLL/info
                scale[j] <- scale[j] - incr
                if (scale[j] <= 0)
                    stop("Parameter estimates (scale) out of bounds")
                if (all(abs(incr) < 1e-05)) 
                    break
            }
            if (iter == maxiter)
                warning("Maximum iterations reached")
        }
        return(list(scale = scale))
    }
    stop("Invalid specification of parameters")
}
