checkbest <- function(y,x,approximate, rhoOpt, resrw, bestscales, besttauscales, c1,c2,b1, bestbetas, Mscale, scalerw, rr, rworst, worsti, betarw) {
    # modified subfunction for the Fast-Tau algorithm for linear regression originally published in 
    # Salibian-Barrera, M., Willems, G. und Zamar, R. (2008): The fast-tau estimator for regression.
    # Journal of Computational and Graphical Statistics, 17 Nr. 3, 659-682.
    if (!approximate) { # compute actual scale, but use tau-conditions!
        scaletest1 <- mean(rhoOpt(resrw / bestscales[worsti],c1)) < b1
        scaletest2 <- sum(rhoOpt(resrw / bestscales[worsti],c2)) < sum(rhoOpt(rworst/bestscales[worsti],c2))
        if (scaletest1 || scaletest2) { # if conditions fulfilled, compute objective value
            snew <- Mscale(resrw, b1, c1, scalerw)
            computeTAU<- TRUE
        }
        else {computeTAU<- FALSE}
    }
    else { # or just compute approximations (and don't bother with the conditions)
        snew <- scalerw
        computeTAU<- TRUE
        if (rr>0) {
            for (kstep in 1:rr) {
                snew <- sqrt( snew^2 * mean( rhoOpt(resrw/snew,c1) ) / b1 )
            }
        }
    }
    if(computeTAU) {
        taunew <- snew * sqrt(mean(rhoOpt(resrw/snew,c2)))
        if (taunew < besttauscales[worsti]) {
            # if candidate has indeed better tau than the worst of the tt best until now, keep it. 
            besttauscales[worsti] <- taunew
            bestscales[worsti] <- snew
            bestbetas[,worsti] <- betarw
            worsti <- which.max(besttauscales)
            rworst <- y - x %*% bestbetas[,worsti]
        }
    }
    results<- list(bestscales=bestscales, besttauscales=besttauscales, bestbetas=bestbetas, rworst=rworst, worsti=worsti)
    return(results)
}
