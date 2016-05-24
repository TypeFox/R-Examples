dw <-
function (x, y = NULL, eps) 
{
    X <- xy.coords(x, y)
    sample <- cbind(X$x, X$y)
    eps.ext <- ashape(sample, alpha = eps)$alpha.extremes
    m <- length(eps.ext)
    n <- dim(sample)[1]
    d <- matrix(0, nrow = m, ncol = n)
    v.x <- matrix(0, nrow = m, ncol = n)
    v.y <- matrix(0, nrow = m, ncol = n)
    theta <- matrix(0, nrow = m, ncol = n)
    for (i in 1:m) {
        v.x[i, ] <- (sample[, 1] - sample[eps.ext[i], 1])
        v.y[i, ] <- (sample[, 2] - sample[eps.ext[i], 2])
        d[i, ] <- sqrt(v.x[i, ]^2 + v.y[i, ]^2)
    }
    v.x[d > 0] <- v.x[d > 0]/d[d > 0]
    v.y[d > 0] <- v.y[d > 0]/d[d > 0]
    theta[d > 0 & d <= 2 * eps] <- acos(d[d > 0 & d <= 2 * eps] * 
        0.5/eps)
    m1 <- cbind(sample[eps.ext, ], eps, 1, 0, pi/2, 1:m, eps.ext, 
        0)
    m2 <- cbind(sample[eps.ext, ], eps, -1, 0, pi/2, 1:m, eps.ext, 
        0)
    archs <- rbind(m1, m2)
    colnames(archs) <- c("c1", "c2", "r", "v.x", "v.y", "theta", 
        "point", "orig.point", "inte")
    n.arc <- dim(archs)[1]
    watch <- 1
    nowatch <- 0
    if (n.arc > 0) {
        while (watch <= n.arc) {
            pwatch <- archs[watch, "point"]
            points <- which(0 < d[pwatch, ] & d[pwatch, ] <= 
                2 * eps)
            if (length(points) > 0) {
                for (j in points) {
                  v.arch <- c(archs[watch, "v.x"], archs[watch, 
                    "v.y"])
                  if (v.arch[2] >= 0) {
                    ang.OX <- acos(v.arch[1])
                  }
                  else {
                    ang.OX <- 2 * pi - acos(v.arch[1])
                  }
                  v.arch.rot <- rotation(v.arch, ang.OX)
                  v.int <- c(v.x[pwatch, j], v.y[pwatch, j])
                  v.int.rot <- rotation(v.int, ang.OX)
                  inter.theta <- theta[pwatch, j]
                  if (v.int.rot[2] >= 0) {
                    ang.v.int.rot.OX <- acos(v.int.rot[1])
                    angles <- c(-archs[watch, "theta"], archs[watch, 
                      "theta"], ang.v.int.rot.OX - inter.theta, 
                      ang.v.int.rot.OX + inter.theta)
                    names(angles) <- c("theta1", "theta2", "beta1", 
                      "beta2")
                    order <- names(sort(angles))
                    theta1 <- ang.OX - archs[watch, "theta"]
                    theta2 <- ang.OX + archs[watch, "theta"]
                    beta1 <- ang.v.int.rot.OX + ang.OX - inter.theta
                    beta2 <- ang.v.int.rot.OX + ang.OX + inter.theta
                    if (all(order == c("theta1", "theta2", "beta1", 
                      "beta2")) | all(order == c("beta1", "beta2", 
                      "theta1", "theta2"))) {
                    }
                    else if (all(order == c("theta1", "beta1", 
                      "theta2", "beta2"))) {
                      ang.middle <- (angles["beta1"] - angles["theta1"])/2
                      v.new <- rotation(c(1, 0), archs[watch, 
                        "theta"] - ang.middle - ang.OX)
                      archs[watch, 1:7] <- c(archs[watch, 1], 
                        archs[watch, 2], archs[watch, 3], v.new[1], 
                        v.new[2], ang.middle, pwatch)
                      archs[watch, 9] <- j
                    }
                    else if (all(order == c("beta1", "theta1", 
                      "theta2", "beta2"))) {
                      nowatch <- c(nowatch, watch)
                    }
                    else if (all(order == c("theta1", "beta1", 
                      "beta2", "theta2"))) {
                      ang.middle <- (angles["beta1"] - angles["theta1"])/2
                      v.new <- rotation(c(1, 0), archs[watch, 
                        "theta"] - ang.middle - ang.OX)
                      ang.middle2 <- (angles["theta2"] - angles["beta2"])/2
                      v.new2 <- rotation(c(1, 0), -archs[watch, 
                        "theta"] + ang.middle2 - ang.OX)
                      archs <- rbind(archs, c(archs[watch, 1], 
                        archs[watch, 2], archs[watch, 3], v.new2[1], 
                        v.new2[2], ang.middle2, pwatch, archs[watch, 
                          8], j))
                      archs[watch, 1:7] <- c(archs[watch, 1], 
                        archs[watch, 2], archs[watch, 3], v.new[1], 
                        v.new[2], ang.middle, pwatch)
                      archs[watch, 9] <- j
                      n.arc <- n.arc + 1
                    }
                    else if (all(order == c("beta1", "theta1", 
                      "beta2", "theta2"))) {
                      ang.middle2 <- (angles["theta2"] - angles["beta2"])/2
                      v.new2 <- rotation(c(1, 0), -archs[watch, 
                        "theta"] + ang.middle2 - ang.OX)
                      archs[watch, 1:7] <- c(archs[watch, 1], 
                        archs[watch, 2], archs[watch, 3], v.new2[1], 
                        v.new2[2], ang.middle2, pwatch)
                      archs[watch, 9] <- j
                    }
                  }
                  else {
                    ang.v.int.rot.OX <- acos(v.int.rot[1])
                    angles <- c(-archs[watch, "theta"], archs[watch, 
                      "theta"], -ang.v.int.rot.OX - inter.theta, 
                      -ang.v.int.rot.OX + inter.theta)
                    names(angles) <- c("theta1", "theta2", "beta1", 
                      "beta2")
                    order <- names(sort(angles))
                    theta1 <- ang.OX - archs[watch, "theta"]
                    theta2 <- ang.OX + archs[watch, "theta"]
                    beta1 <- -ang.v.int.rot.OX + ang.OX - inter.theta
                    beta2 <- -ang.v.int.rot.OX + ang.OX + inter.theta
                    if (all(order == c("theta1", "theta2", "beta1", 
                      "beta2")) | all(order == c("beta1", "beta2", 
                      "theta1", "theta2"))) {
                    }
                    else if (all(order == c("theta1", "beta1", 
                      "theta2", "beta2"))) {
                      ang.middle <- (angles["beta1"] - angles["theta1"])/2
                      v.new <- rotation(c(1, 0), archs[watch, 
                        "theta"] - ang.middle - ang.OX)
                      archs[watch, 1:7] <- c(archs[watch, 1], 
                        archs[watch, 2], archs[watch, 3], v.new[1], 
                        v.new[2], ang.middle, pwatch)
                      archs[watch, 9] <- j
                    }
                    else if (all(order == c("beta1", "theta1", 
                      "theta2", "beta2"))) {
                      nowatch <- c(nowatch, watch)
                    }
                    else if (all(order == c("theta1", "beta1", 
                      "beta2", "theta2"))) {
                      ang.middle <- (angles["beta1"] - angles["theta1"])/2
                      v.new <- rotation(c(1, 0), archs[watch, 
                        "theta"] - ang.middle - ang.OX)
                      ang.middle2 <- (angles["theta2"] - angles["beta2"])/2
                      v.new2 <- rotation(c(1, 0), -archs[watch, 
                        "theta"] + ang.middle2 - ang.OX)
                      archs <- rbind(archs, c(archs[watch, 1], 
                        archs[watch, 2], archs[watch, 3], v.new2[1], 
                        v.new2[2], ang.middle2, pwatch, archs[watch, 
                          8], j))
                      archs[watch, 1:7] <- c(archs[watch, 1], 
                        archs[watch, 2], archs[watch, 3], v.new[1], 
                        v.new[2], ang.middle, pwatch)
                      archs[watch, 9] <- j
                      n.arc <- n.arc + 1
                    }
                    else if (all(order == c("beta1", "theta1", 
                      "beta2", "theta2"))) {
                      ang.middle2 <- (angles["theta2"] - angles["beta2"])/2
                      v.new2 <- rotation(c(1, 0), -archs[watch, 
                        "theta"] + ang.middle2 - ang.OX)
                      archs[watch, 1:7] <- c(archs[watch, 1], 
                        archs[watch, 2], archs[watch, 3], v.new2[1], 
                        v.new2[2], ang.middle2, pwatch)
                      archs[watch, 9] <- j
                    }
                  }
                }
            }
            watch <- watch + 1
        }
    }
    if (length(nowatch) > 1) {
        archs <- archs[-nowatch, ]
    }
    return(archs)
}
