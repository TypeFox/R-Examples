gammaFit <- function(X, x, y, maxiter=100, show.iter=FALSE, tol=0.001){
    g.old <- 1
    converged <- FALSE
    step <- 1
    v <- rep(1, length(x))

    while(!converged & (step <= maxiter)) {
        g.hat <- gamEst(X, x, y, v)
        if (show.iter) cat("step",step,"gamma = ",g.hat,"\n")
        v <- x^g.hat
        if(abs((g.old - g.hat)/g.old) < tol) {converged <- TRUE}
        g.old <- g.hat
        v <- x^g.hat
        step <- step+1
    }

    if ((step >= maxiter) & !converged){
        cat("Maximum no. of iterations reached without convergence.\n")
        cat("g.hat = ", g.hat, "\n")
    }
    else{
        cat("Convergence attained in ", step-1, "steps.\n")
        cat("g.hat =", g.hat, "\n")
    }
    list(g.hat=g.hat, converged = converged, steps=step-1)
}