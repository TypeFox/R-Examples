monmlp.initialize <-
function(x, y, hidden1, hidden2, init.weights)
{
    w.min <- min(init.weights)
    w.max <- max(init.weights)
    W1 <- matrix(runif((ncol(x)+1)*hidden1, w.min, w.max),
                 ncol(x)+1, hidden1)
    if (hidden2==0){
        W2 <- matrix(runif((hidden1+1)*ncol(y), w.min, w.max),
                     hidden1+1, ncol(y))
        W.vec <- c(W1, W2)
    } else{
        W2 <- matrix(runif((hidden1+1)*hidden2, w.min, w.max),
                     hidden1+1, hidden2)
        W3 <- matrix(runif((hidden2+1)*ncol(y), w.min, w.max),
                     hidden2+1, ncol(y))
        W.vec <- c(W1, W2, W3)
    }
    W.vec
}

