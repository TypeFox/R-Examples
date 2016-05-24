# Set method to predict Author: Yiwen Zhang

setMethod("predict", "MGLMreg", function(object, newdata) {
    beta <- object$coefficients
    d <- ncol(object$data$Y)
    dist <- object$distribution
    
    if (dist == "Multinomial") {
        pred <- exp(newdata %*% beta)/(rowSums(exp(newdata %*% beta)) + 1)
        pred <- cbind(pred, (1 - rowSums(pred)))
    } else if (dist == "Dirichlet Multinomial") {
        pred <- exp(newdata %*% beta)/rowSums(exp(newdata %*% beta))
    } else if (dist == "Generalized Dirichlet Multinomial") {
        alpha <- beta[, 1:(d - 1)]
        beta <- beta[, d:ncol(beta)]
        pred <- matrix(0, nrow(newdata), d)
        pred[, 1:(d - 1)] <- exp(newdata %*% alpha)/(exp(newdata %*% alpha) + exp(newdata %*% 
            beta))
        pred[, 2] <- (1 - pred[, 1]) * pred[, 2]
        if (nrow(newdata) > 1) {
            for (f in 3:(d - 1)) {
                pred[, f] <- (1 - rowSums(pred[, 1:(f - 1)])) * pred[, f]
            }
            pred[, d] <- (1 - rowSums(pred[, 1:(d - 1)]))
        } else {
            for (f in 3:(d - 1)) {
                pred[, f] <- (1 - sum(pred[, 1:(f - 1)])) * pred[, f]
            }
            pred[, d] <- (1 - sum(pred[, 1:(d - 1)]))
        }
    } else if (dist == "Negative Multinomial") {
        if (is.null(beta$phi)) {
            alpha <- beta[, 1:d]
            beta <- beta[, (d + 1)]
        } else {
            alpha <- beta$alpha
            phi <- beta$phi
        }
        pred <- exp(newdata %*% alpha)/(rowSums(exp(newdata %*% alpha)) + 1)
    }
    colnames(pred) <- colnames(object$data$Y)
    return(pred)
}) 
