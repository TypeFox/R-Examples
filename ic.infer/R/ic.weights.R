ic.weights <- function (corr, ...) 
{
    ## ... given to pmvnorm, should contain algorithm information only
    ## changed to ... in order to stay tuned with updates in mvtnorm
    ## and be also compatible with old versions (e.g. 2.5.0)
    if (!(is.matrix(corr) & nrow(corr) == ncol(corr))) 
        stop("corr must be a square matrix.")
    if (!(all(eigen(corr)$value > 0))) 
        stop("corr must be positive definite.")
    g <- nrow(corr)
    liste <- 1:g
    ## weights in the order of highest to lowest dimension 
    ## for variable with correlation or cov corr
    weights <- rep(0, g + 1)
    names(weights) <- rev(0:g)
    weights[1] <- pmvnorm(lower = rep(0, g), upper = rep(Inf, 
        g), sigma = corr, ...)
    weights[g + 1] <- pmvnorm(lower = rep(0, g), upper = rep(Inf, 
        g), sigma = solve(corr), ...)
    ## prevent lengthy calculations 
    ## if longest vector too long for available storage
    if (g > 4){
    if (!is.numeric(try(matrix(0, floor((g-2)/2), choose(g, floor((g-2)/2))), 
        silent = TRUE))) 
        stop(paste("ic.weights will not work, corr too large, \n", 
            "interim matrix with ", floor((g-2)/2) * choose(g, floor((g-2)/2)), 
            " elements cannot be created.", sep = ""))
    }
    if (g==2) weights[2] <- 1-sum(weights)
    if (g==3) {
        weights[2] <- 0.5 - weights[4]
        weights[3] <- 0.5 - weights[1]
    }
    if (g > 3) {
        for (k in 1:floor((g - 2)/2)) {
            ### fill weights simultaneously from top and bottom
            ### fill the two middle weights by 0.5 - odd and even sum of others, 
            ###     respectively, according to Silvapulle and Sen 2004, 
            ###     Prop. 3.6.1, part 3
            jetzt <- nchoosek(g, k)
            wjetzt <- matrix(0, choose(g, k), 2)
            for (j in 1:(choose(g, k))) {
                diese <- jetzt[, j]
                andere <- setdiff(liste, diese)
                hilf <- corr[andere, andere] - corr[andere, diese] %*% 
                  solve(corr[diese, diese], matrix(corr[diese, 
                    andere], k, g - k))
                    ## matrix necessary, because vector otherwise wrong direction
                wjetzt[j, 1] <- pmvnorm(lower = rep(0, k), upper = rep(Inf, 
                  k), sigma = solve(corr[diese, diese]), ...) * 
                  pmvnorm(lower = rep(0, g - k), upper = rep(Inf, 
                    g - k), sigma = hilf, ...)
                if (!k==(g-2)/2){
                  ## needed for odd g only for the last one from top, 
                  ## for even g calculated as 0.5 minus the other even weights
                  hilf <- corr[diese, diese] - corr[diese, andere] %*% 
                    solve(corr[andere, andere], matrix(corr[andere, 
                      diese], g - k, k))
                  wjetzt[j, 2] <- pmvnorm(lower = rep(0, g - k), 
                    upper = rep(Inf, g - k), sigma = solve(corr[andere, 
                      andere]), ...) * pmvnorm(lower = rep(0, k), 
                    upper = rep(Inf, k), sigma = hilf, ...)
                }
            }
            weights[k + 1] <- sum(wjetzt[, 1])
            weights[g + 1 - k] <- sum(wjetzt[, 2])
        }
    ### fill last weight using sum of odd and even weights 
    ### even weights in odd positions and vice versa
    if (g/2 == floor(g/2)) {
        even.sum <- sum(weights[1+2*((g/2):0)])
        odd.sum <- sum(weights[2*((g/2):1)])
        if (g/4 == floor(g/4)) {
          weights[g/2 + 1] <- 0.5 - even.sum ## even weights
          weights[g/2 + 2] <- 0.5 - odd.sum ## odd weights
        }
        else {
          weights[g/2 + 1] <- 0.5 - odd.sum ## odd weights
          weights[g/2 + 2] <- 0.5 - even.sum ## even weights
        }
    }
    else {
        even.sum <- sum(weights[2*(((g+1)/2):1)])
        odd.sum <- sum(weights[2*(((g+1)/2):1)-1])
        if ((g+1)/4 == floor((g+1)/4)) {
          weights[(g + 1)/2] <- 0.5 - even.sum ## even weights
          weights[(g + 3)/2] <- 0.5 - odd.sum ## odd weights
        }
        else {
          weights[(g + 1)/2] <- 0.5 - odd.sum ## odd weights
          weights[(g + 3)/2] <- 0.5 - even.sum ## even weights
        }
    }
    }
    weights
}
