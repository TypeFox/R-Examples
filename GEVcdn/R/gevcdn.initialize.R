gevcdn.initialize <-
function (x, n.hidden, init.range)
{
    W1 <- matrix(runif((ncol(x) + 1)*n.hidden, min = min(init.range),
                 max = max(init.range)), nrow = ncol(x) + 1,
                 ncol = n.hidden)
    W2 <- matrix(runif((n.hidden + 1)*3, min = min(init.range),
                 max = max(init.range)), nrow = n.hidden + 1, ncol = 3)
    c(W1, W2)
}

