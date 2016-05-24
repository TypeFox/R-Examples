confidence.interval <-
function (x, alpha = 0.05) 
{
    net <- x
    covariate <- cbind(1, net$covariate)
    response <- net$response
    err.fct <- net$err.fct
    act.fct <- net$act.fct
    linear.output <- net$linear.output
    exclude <- net$exclude
    list.weights <- net$weights
    rep <- length(list.weights)
    if (is.null(list.weights)) 
        stop("weights were not calculated")
    nrow.weights <- sapply(list.weights[[1]], nrow)
    ncol.weights <- sapply(list.weights[[1]], ncol)
    lower.ci <- NULL
    upper.ci <- NULL
    nic <- NULL
    for (i in 1:rep) {
        weights <- list.weights[[i]]
        error <- net$result.matrix["error", i]
        if (length(weights) > 2) 
            stop("nic and confidence intervals will not be calculated for more than one hidden layer of neurons", 
                call. = FALSE)
        result.nic <- calculate.information.matrices(covariate, 
            response, weights, err.fct, act.fct, exclude, linear.output)
        nic <- c(nic, error + result.nic$trace)
        if (!is.null(result.nic$variance)) {
            if (all(diag(result.nic$variance) >= 0)) {
                weights.vector <- unlist(weights)
                if (!is.null(exclude)) {
                  d <- rep(NA, length(weights.vector))
                  d[-exclude] <- qnorm(1 - alpha/2) * sqrt(diag(result.nic$variance))/sqrt(nrow(covariate))
                }
                else {
                  d <- qnorm(1 - alpha/2) * sqrt(diag(result.nic$variance))/sqrt(nrow(covariate))
                }
                lower.ci <- c(lower.ci, list(relist(weights.vector - 
                  d, nrow.weights, ncol.weights)))
                upper.ci <- c(upper.ci, list(relist(weights.vector + 
                  d, nrow.weights, ncol.weights)))
            }
        }
    }
    if (length(lower.ci) < rep) 
        warning(sprintf("%s of %s repetition(s) could not calculate confidence intervals for the weights; varify that the neural network does not contain irrelevant neurons", 
            length(lower.ci), rep), call. = F)
    list(lower.ci = lower.ci, upper.ci = upper.ci, nic = nic)
}
calculate.information.matrices <-
function (covariate, response, weights, err.fct, act.fct, exclude, 
    linear.output) 
{
    temp <- act.fct
    if (type(act.fct) == "logistic") {
        act.deriv.fct <- function(x) {
            act.fct(x) * (1 - act.fct(x))
        }
        act.deriv2.fct <- function(x) {
            act.fct(x) * (1 - act.fct(x)) * (1 - (2 * act.fct(x)))
        }
    }
    else {
        attr(temp, "type") <- NULL
        act.deriv.fct <- differentiate(temp)
        act.deriv2.fct <- differentiate(temp, hessian = T)
    }
    temp <- err.fct
    attr(temp, "type") <- NULL
    err.deriv.fct <- differentiate(temp)
    err.deriv2.fct <- differentiate(temp, hessian = T)
    length.weights <- length(weights)
    nrow.weights <- sapply(weights, nrow)
    ncol.weights <- sapply(weights, ncol)
    if (linear.output) {
        output.act.fct <- function(x) {
            x
        }
        output.act.deriv.fct <- function(x) {
            matrix(1, nrow(x), ncol(x))
        }
        output.act.deriv2.fct <- function(x) {
            matrix(0, nrow(x), ncol(x))
        }
    }
    else {
        output.act.fct <- act.fct
        output.act.deriv.fct <- act.deriv.fct
        output.act.deriv2.fct <- act.deriv2.fct
    }
    neuron.deriv <- NULL
    neuron.deriv2 <- NULL
    neurons <- list(covariate)
    if (length.weights > 1) 
        for (i in 1:(length.weights - 1)) {
            temp <- neurons[[i]] %*% weights[[i]]
            act.temp <- act.fct(temp)
            neuron.deriv[[i]] <- act.deriv.fct(temp)
            neuron.deriv2[[i]] <- act.deriv2.fct(temp)
            neurons[[i + 1]] <- cbind(1, act.temp)
        }
    if (!is.list(neuron.deriv)) 
        neuron.deriv <- list(neuron.deriv)
    if (!is.list(neuron.deriv2)) 
        neuron.deriv2 <- list(neuron.deriv2)
    temp <- neurons[[length.weights]] %*% weights[[length.weights]]
    net.result <- output.act.fct(temp)
    neuron.deriv[[length.weights]] <- output.act.deriv.fct(temp)
    neuron.deriv2[[length.weights]] <- output.act.deriv2.fct(temp)
    err.deriv <- err.deriv.fct(net.result, response)
    err.deriv2 <- err.deriv2.fct(net.result, response)
    if (any(is.na(unlist(neuron.deriv2)))) {
        return(list(nic = NA, hessian = NULL))
        warning("neuron.deriv2 contains NA; this might be caused by a wrong choice of 'act.fct'")
    }
    if (any(is.na(err.deriv)) || any(is.na(err.deriv2))) {
        if (type(err.fct) == "ce") {
            one <- which(net.result == 1)
            if (length(one) > 0) {
                for (i in 1:length(one)) {
                  if (response[one[i]] == 1) {
                    err.deriv[one[i]] <- 1
                    err.deriv2[one[i]] <- 1
                  }
                }
            }
            zero <- which(net.result == 0)
            if (length(zero) > 0) {
                for (i in 1:length(zero)) {
                  if (response[zero[i]] == 0) {
                    err.deriv[zero[i]] <- 1
                    err.deriv2[zero[i]] <- -1
                  }
                }
            }
        }
    }
    if (any(is.na(err.deriv))) {
        return(list(nic = NA, hessian = NULL))
        warning("err.deriv contains NA; this might be caused by a wrong choice of 'act.fct'")
    }
    if (any(is.na(err.deriv2))) {
        return(list(nic = NA, hessian = NULL))
        warning("err.deriv2 contains NA; this might be caused by a wrong choice of 'act.fct'")
    }
    if (length.weights == 2) {
        length.betha <- (nrow.weights * ncol.weights)[1]
        length.alpha <- (nrow.weights * ncol.weights)[2]
        total.length.weights <- length.alpha + length.betha
        betha.ind <- matrix(1:length.betha, nrow = nrow.weights[1], 
            ncol = ncol.weights[1])
        alpha.ind <- matrix(1:length.alpha, nrow = nrow.weights[2], 
            ncol = ncol.weights[2])
        Hesse <- matrix(NA, nrow = total.length.weights, ncol = total.length.weights)
        Cross.Gradient <- matrix(NA, nrow = total.length.weights, 
            ncol = total.length.weights)
        Cross.Gradient2 <- matrix(NA, nrow = total.length.weights, 
            ncol = total.length.weights)
        for (i in 1:total.length.weights) {
            for (j in 1:total.length.weights) {
                if (is.null(exclude) || all(i != exclude & j != 
                  exclude)) {
                  if (i <= length.betha) {
                    temp <- which(betha.ind == i, arr.ind = T)
                    r <- temp[1]
                    s <- temp[2]
                  }
                  else {
                    temp <- which(alpha.ind == (i - length.betha), 
                      arr.ind = T)
                    r <- temp[1]
                    s <- temp[2]
                  }
                  if (j <= length.betha) {
                    temp <- which(betha.ind == j, arr.ind = T)
                    u <- temp[1]
                    v <- temp[2]
                  }
                  else {
                    temp <- which(alpha.ind == (j - length.betha), 
                      arr.ind = T)
                    u <- temp[1]
                    v <- temp[2]
                  }
                  if ((i <= length.betha) && (j <= length.betha)) {
                    Cross.Gradient[i, j] <- sum(((err.deriv^2 * 
                      neuron.deriv[[2]]^2) %*% (weights[[2]][(s + 
                      1), ] * weights[[2]][(v + 1), ])) * neuron.deriv[[1]][, 
                      s] * neurons[[1]][, r] * neuron.deriv[[1]][, 
                      v] * neurons[[1]][, u])
                    Cross.Gradient2[i, j] <- sum(((err.deriv2 * 
                      neuron.deriv[[2]]^2) %*% (weights[[2]][(s + 
                      1), ] * weights[[2]][(v + 1), ])) * neuron.deriv[[1]][, 
                      s] * neurons[[1]][, r] * neuron.deriv[[1]][, 
                      v] * neurons[[1]][, u])
                    if (s == v) 
                      Hesse[i, j] <- sum(neuron.deriv[[1]][, 
                        s] * neurons[[1]][, r] * neuron.deriv[[1]][, 
                        v] * neurons[[1]][, u] * ((neuron.deriv2[[2]] * 
                        err.deriv) %*% (weights[[2]][(s + 1), 
                        ] * weights[[2]][(v + 1), ]))) + sum(neuron.deriv2[[1]][, 
                        s] * neurons[[1]][, r] * neurons[[1]][, 
                        u] * ((neuron.deriv[[2]] * err.deriv) %*% 
                        weights[[2]][(s + 1), ]))
                    else Hesse[i, j] <- sum(neuron.deriv[[1]][, 
                      s] * neurons[[1]][, r] * neuron.deriv[[1]][, 
                      v] * neurons[[1]][, u] * ((neuron.deriv2[[2]] * 
                      err.deriv) %*% (weights[[2]][(s + 1), ] * 
                      weights[[2]][(v + 1), ])))
                  }
                  else if ((i > length.betha) && (j > length.betha)) {
                    if (v == s) {
                      Cross.Gradient[i, j] <- sum(err.deriv[, 
                        v]^2 * (neuron.deriv[[2]][, s] * neurons[[2]][, 
                        r] * neuron.deriv[[2]][, v] * neurons[[2]][, 
                        u]))
                      Cross.Gradient2[i, j] <- sum(err.deriv2[, 
                        v] * (neuron.deriv[[2]][, s] * neurons[[2]][, 
                        r] * neuron.deriv[[2]][, v] * neurons[[2]][, 
                        u]))
                    }
                    else {
                      Cross.Gradient[i, j] <- 0
                      Cross.Gradient2[i, j] <- 0
                    }
                    if (v == s) 
                      Hesse[i, j] <- sum(neuron.deriv2[[2]][, 
                        s] * err.deriv[, s] * neurons[[2]][, 
                        u] * neurons[[2]][, r])
                    else Hesse[i, j] <- 0
                  }
                  else if ((i > length.betha) && (j <= length.betha)) {
                    Cross.Gradient[i, j] <- sum(err.deriv[, s]^2 * 
                      (neuron.deriv[[2]][, s] * neurons[[2]][, 
                        r] * (neuron.deriv[[2]][, s] * weights[[2]][(v + 
                        1), s]) * neuron.deriv[[1]][, v] * neurons[[1]][, 
                        u]))
                    Cross.Gradient2[i, j] <- sum(err.deriv2[, 
                      s] * (neuron.deriv[[2]][, s] * neurons[[2]][, 
                      r] * (neuron.deriv[[2]][, s] * weights[[2]][(v + 
                      1), s]) * neuron.deriv[[1]][, v] * neurons[[1]][, 
                      u]))
                    if (v == r) 
                      Hesse[i, j] <- sum(neurons[[2]][, r] * 
                        neuron.deriv[[1]][, v] * neurons[[1]][, 
                        u] * neuron.deriv2[[2]][, s] * err.deriv[, 
                        s] * weights[[2]][(v + 1), s]) + sum(neuron.deriv[[2]][, 
                        s] * err.deriv[, s] * neurons[[1]][, 
                        u] * neuron.deriv[[1]][, v])
                    else Hesse[i, j] <- sum(neurons[[2]][, r] * 
                      neuron.deriv[[1]][, v] * neurons[[1]][, 
                      u] * neuron.deriv2[[2]][, s] * err.deriv[, 
                      s] * weights[[2]][(v + 1), s])
                  }
                  else {
                    Cross.Gradient[i, j] <- sum(err.deriv[, v]^2 * 
                      (neuron.deriv[[2]][, v] * neurons[[2]][, 
                        u] * (neuron.deriv[[2]][, v] * weights[[2]][(s + 
                        1), v]) * neuron.deriv[[1]][, s] * neurons[[1]][, 
                        r]))
                    Cross.Gradient2[i, j] <- sum(err.deriv2[, 
                      v] * (neuron.deriv[[2]][, v] * neurons[[2]][, 
                      u] * (neuron.deriv[[2]][, v] * weights[[2]][(s + 
                      1), v]) * neuron.deriv[[1]][, s] * neurons[[1]][, 
                      r]))
                    if (s == u) 
                      Hesse[i, j] <- sum(neurons[[2]][, u] * 
                        neuron.deriv[[1]][, s] * neurons[[1]][, 
                        r] * neuron.deriv2[[2]][, v] * err.deriv[, 
                        v] * weights[[2]][(s + 1), v]) + sum(neuron.deriv[[2]][, 
                        v] * err.deriv[, v] * neurons[[1]][, 
                        r] * neuron.deriv[[1]][, s])
                    else Hesse[i, j] <- sum(neurons[[2]][, u] * 
                      neuron.deriv[[1]][, s] * neurons[[1]][, 
                      r] * neuron.deriv2[[2]][, v] * err.deriv[, 
                      v] * weights[[2]][(s + 1), v])
                  }
                }
            }
        }
    }
    else if (length.weights == 1) {
        length.alpha <- sum(nrow.weights * ncol.weights)
        alpha.ind <- matrix(1:length.alpha, nrow = nrow.weights[1], 
            ncol = ncol.weights[1])
        Hesse <- matrix(NA, nrow = length.alpha, ncol = length.alpha)
        Cross.Gradient <- matrix(NA, nrow = length.alpha, ncol = length.alpha)
        Cross.Gradient2 <- matrix(NA, nrow = length.alpha, ncol = length.alpha)
        for (i in 1:length.alpha) {
            for (j in 1:length.alpha) {
                if (is.null(exclude) || all(i != exclude & j != 
                  exclude)) {
                  r <- which(alpha.ind == i, arr.ind = T)[1]
                  s <- which(alpha.ind == i, arr.ind = T)[2]
                  u <- which(alpha.ind == j, arr.ind = T)[1]
                  v <- which(alpha.ind == j, arr.ind = T)[2]
                  if (s == v) {
                    Hesse[i, j] <- sum(neuron.deriv2[[1]][, s] * 
                      err.deriv[, s] * neurons[[1]][, r] * neurons[[1]][, 
                      u])
                    Cross.Gradient[i, j] <- sum(neuron.deriv[[1]][, 
                      s]^2 * err.deriv[, s]^2 * neurons[[1]][, 
                      r] * neurons[[1]][, u])
                    Cross.Gradient2[i, j] <- sum(neuron.deriv[[1]][, 
                      s]^2 * err.deriv2[, s] * neurons[[1]][, 
                      r] * neurons[[1]][, u])
                  }
                  else {
                    Hesse[i, j] <- 0
                    Cross.Gradient[i, j] <- 0
                    Cross.Gradient2[i, j] <- 0
                  }
                }
            }
        }
    }
    B <- Cross.Gradient/nrow(neurons[[1]])
    A <- (Cross.Gradient2 + Hesse)/nrow(neurons[[1]])
    if (!is.null(exclude)) {
        B <- as.matrix(B[-exclude, -exclude])
        A <- as.matrix(A[-exclude, -exclude])
    }
    if (det(A) == 0) {
        trace <- NA
        variance <- NULL
    }
    else {
        A.inv <- ginv(A)
        variance <- A.inv %*% B %*% A.inv
        trace <- sum(diag(B %*% A.inv))
    }
    return(list(trace = trace, variance = variance))
}
