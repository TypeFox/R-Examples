compute <-
function (x, covariate, rep = 1) 
{
    nn <- x
    linear.output <- nn$linear.output
    weights <- nn$weights[[rep]]
    nrow.weights <- sapply(weights, nrow)
    ncol.weights <- sapply(weights, ncol)
    weights <- unlist(weights)
    if (any(is.na(weights))) 
        weights[is.na(weights)] <- 0
    weights <- relist(weights, nrow.weights, ncol.weights)
    length.weights <- length(weights)
    covariate <- as.matrix(cbind(1, covariate))
    act.fct <- nn$act.fct
    neurons <- list(covariate)
    if (length.weights > 1) 
        for (i in 1:(length.weights - 1)) {
            temp <- neurons[[i]] %*% weights[[i]]
            act.temp <- act.fct(temp)
            neurons[[i + 1]] <- cbind(1, act.temp)
        }
    temp <- neurons[[length.weights]] %*% weights[[length.weights]]
    if (linear.output) 
        net.result <- temp
    else net.result <- act.fct(temp)
    list(neurons = neurons, net.result = net.result)
}
