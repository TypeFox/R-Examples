qrnn.eval <-
function(x, W1, W2, lower, eps, Th)
{
    x <- cbind(x, 1)
    h1 <- x %*% W1
    y1 <- Th(h1)
    aug.y1 <- cbind(y1, 1)
    y2 <- aug.y1 %*% W2
    y2 <- hramp(y2, lower, eps)
    y2
}

