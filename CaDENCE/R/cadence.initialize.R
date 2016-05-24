cadence.initialize <-
function(x, n.hidden, init.range, distribution)
{
    n.parms <- length(distribution$parameters)
    W1 <- matrix(runif((ncol(x)+1)*n.hidden, init.range[1], init.range[2]),
                 ncol(x)+1, n.hidden)
    W2 <- matrix(runif((n.hidden+1)*n.parms, init.range[1], init.range[2]),
                 n.hidden+1, n.parms)
    c(W1, W2)
}

