Sample <-
function(channel, n, index, show.progress = FALSE) {
    ## Generate random samples according to the distribution P(y | x)
    ## described by a given channel.
    
    x.out <- matrix(data = NA, nrow = n, ncol = ncol(channel$x))
    y.out <- matrix(data = NA, nrow = n, ncol = ncol(channel$y))

    if(show.progress) {
        bar <- txtProgressBar(min = 0, max = n,
                              width = 40, style = 3, char = "*")
    }
    for(i in 1:n) {
        if(show.progress) setTxtProgressBar(bar, i)
        
        x.out[i, ] <- channel$x[index, ]

        cpd <- ConditionalDistribution(channel, index)
        k <- which(rmultinom(1, 1, cpd$p) == 1)
        y.out[i, ] <- cpd$y[k, ]
    }
    if(show.progress) cat("\n")
    
    list(x = x.out, y = y.out)
}
