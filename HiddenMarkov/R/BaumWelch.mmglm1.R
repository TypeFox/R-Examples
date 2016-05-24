BaumWelch.mmglm1 <- function (object, control=bwcontrol(), ...){
    tol <- control$tol
    oldLL <- -Inf
    for (iter in 1:control$maxiter) {
        cond <- Estep.mmglm1(object)
        diff <- cond$LL - oldLL
        if (control$prt) {
            cat("iter =", iter, "\n")
            cat("LL =", formatC(cond$LL, digits=log10(1/tol)+2,
                                format="f"), "\n")
            cat("diff =", diff, "\n\n")
        }
        if (diff < 0 & control$posdiff)
            stop("Worse log-likelihood on last iteration")
        if (eval(control$converge)) break
        #----  Mstep  ----
        Pi <- diag(1/apply(cond$v, MARGIN=2, FUN=sum)) %*% 
                  apply(cond$v, MARGIN=c(2, 3), FUN=sum)
        if (object$nonstat) delta <- cond$u[1, ]
            else            delta <- compdelta(Pi)

        tmp <- Mstep.mmglm1(object, cond$u)
        #-----------------
        oldLL <- cond$LL
        object$delta <- delta
        object$Pi <- Pi
        object$beta <- tmp$beta
        object$sigma <- tmp$sigma
    }
    rownames(object$beta) <- colnames(object$Xdesign)
    colnames(object$beta) <- paste("State", 1:length(object$delta))
    object$u <- cond$u    
    object$v <- cond$v
    object$LL <- cond$LL
    object$iter <- iter
    object$diff <- diff
    return(object)
}

