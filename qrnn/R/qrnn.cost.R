qrnn.cost <-
function(weights, x, y, n.hidden, tau, lower, eps, Th, Th.prime,
         penalty)
{
    penalty2 <- ifelse(identical(Th, linear), penalty, 0)
    w <- qrnn.reshape(x, y, weights, n.hidden)
    W1 <- w$W1; rW1 <- nrow(W1); cW1 <- ncol(W1)
    W2 <- w$W2; rW2 <- nrow(W2); cW2 <- ncol(W2)
    # Forward pass
    x <- cbind(x, 1)
    h1 <- x %*% W1
    y1 <- Th(h1)
    aug.y1 <- cbind(y1, 1)
    h2 <- aug.y1 %*% W2
    y2 <- hramp(h2, lower, eps)
    E <- y-y2
    # Backward pass
    delta2 <- hramp.prime(h2, lower, eps)*
              tilted.huber.prime(E, tau, eps)
    gradient.W2 <- -(t(aug.y1) %*% delta2)/length(E) +
        2*penalty2*rbind(W2[1:(rW2-1),,drop=FALSE], 0)/(length(W2)-cW2)
    E1 <- delta2 %*% t(W2[1:(rW2-1),,drop=FALSE])
    delta1 <- Th.prime(h1)*E1
    gradient.W1 = -(t(x) %*% delta1)/length(E) +
        2*penalty*rbind(W1[1:(rW1-1),,drop=FALSE], 0)/(length(W1)-cW1)
    # Error & gradient
    cost <- sum(tilted.huber(E, tau, eps))/length(E) +
         penalty*sum(W1[1:(rW1-1),,drop=FALSE]^2)/(length(W1)-cW1) +
        penalty2*sum(W2[1:(rW2-1),,drop=FALSE]^2)/(length(W2)-cW2)
    gradient <- c(gradient.W1, gradient.W2)
    attr(cost, "gradient") <- gradient
    cost
}

