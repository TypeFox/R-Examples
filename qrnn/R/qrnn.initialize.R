qrnn.initialize <-
function(x, y, n.hidden)
{
    W1 <- matrix(runif((ncol(x)+1)*n.hidden, -0.5, 0.5),
                 ncol(x)+1, n.hidden)
    W2 <- matrix(runif((n.hidden+1)*ncol(y), -0.5, 0.5),
                 n.hidden+1, ncol(y))
    c(W1, W2)
}

