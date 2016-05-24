BaumWelch.mmpp <- function (object, control = bwcontrol(), ...){
    #   scaled version of Ryden (1996)
    tau <- object$tau[-1] - object$tau[-length(object$tau)]
    Q <- object$Q
    delta <- object$delta
    lambda <- object$lambda
    tol <- control$tol
    #---------------------------------------------
    m <- nrow(Q)
    n <- length(tau)
    oldLL <- -Inf
    for (iter in 1:control$maxiter) {
        cond <- Estep.mmpp(tau, Q, delta, lambda)
        LL <- cond$LL
        diff <- LL - oldLL
        if (control$prt) {
            cat("iter =", iter, "\n")
            cat("LL =", formatC(LL, digits=log10(1/tol)+2,
                                format="f"), "\n")
            cat("diff =", diff, "\n\n")
        }
        if (diff < 0 & control$posdiff) stop("Worse log-likelihood on last iteration")
        if (eval(control$converge)) break
        #----  Mstep  ----
        Q <- Q * (diag(1/diag(cond$A)) %*% cond$A)
        diag(Q) <- 0
        diag(Q) <- -apply(Q, MARGIN=1, FUN=sum)
        lambda <- cond$B/diag(cond$A)
        if (object$nonstat) delta <- exp(cond$logalpha[1, ] +
                                  cond$logbeta[1, ])
        else delta <- compdelta(solve(diag(lambda) -
                                Q) %*% diag(lambda))
        oldLL <- LL
    }
    x <- list(delta = delta, Q = Q, lambda = lambda,
              LL = LL, iter = iter, diff = diff)
    #-----------------------------------------------
    nms <- names(x)
    k <- length(nms)
    for (i in 1:k) object[[nms[i]]] <- x[[nms[i]]]
    return(object)
}

