rprop <-
function(w, f, iterlim=100, print.level=1, delta.0=0.1, delta.min=1E-6,
         delta.max=50, epsilon=1e-8, step.tol=1e-6, f.target=-Inf, ...)
{
    fprime <-
    function(w.init, f, epsilon, f.init, ...)
    {
        if(is.null(f.init)) f.init <- f(w.init, ...)
        gradient <- w.init*0
        for(i in 1:length(w.init)){
            w.plus <- w.init
            w.plus[i] <- w.plus[i] + epsilon
            f.plus <- f(w.plus, ...)
            gradient[i] <- (f.plus-f.init)/epsilon
        }
        gradient
    }
    dEdw.prev <- delta.w <- w*0
    delta <- w*0 + delta.0
    nu.minus <- 0.5
    nu.plus <- 1.2
    E.iter <- rep(Inf, iterlim)
    dE <- Inf
    for (i in 1:iterlim){
        E <- f(w, ...)
        dEdw <- attr(E, "gradient")
        if(is.null(dEdw)){
            dEdw <- fprime(w.init=w, f=f, epsilon=epsilon, f.init=E, ...)
        }
        E.iter[i] <- E
        if(i > 3) dE <- mean(abs(diff(E.iter[(i-3):i])))
        if (print.level > 0) cat(i, E, dE, "\n")
        if((dE < step.tol) | (E < f.target)){
            return(list(par=w, value=E, gradient=dEdw))
        }
        i.pos <- (dEdw.prev*dEdw) > 0
        i.neg <- (dEdw.prev*dEdw) < 0
        i.zero <- !(i.pos | i.neg)
        if(any(i.pos)){
            delta[i.pos] <- pmin(delta[i.pos]*nu.plus, delta.max)
            delta.w[i.pos] <- (-sign(dEdw)*delta)[i.pos]
            w[i.pos] <- w[i.pos] + delta.w[i.pos]
            dEdw.prev[i.pos] <- dEdw[i.pos]
        }
        if(any(i.neg)){
            delta[i.neg] <- pmax(delta[i.neg]*nu.minus, delta.min)
            if((i > 1) & (E.iter[i] > E.iter[i-1]))
                w[i.neg] <- w[i.neg] - delta.w[i.neg]
            dEdw.prev[i.neg] <- 0
        }
        if(any(i.zero)){
            delta.w[i.zero] <- (-sign(dEdw)*delta)[i.zero]
            w[i.zero] <- w[i.zero] + delta.w[i.zero]
            dEdw.prev[i.zero] <- dEdw[i.zero]
        }
    }
    E <- f(w, ...)
    dEdw <- attr(E, "gradient")
    if(is.null(dEdw)){
        dEdw <- fprime(w.init=w, f=f, epsilon=epsilon, f.init=E, ...)
    }
    if (print.level > 0) cat("*", E, "\n")
    list(par=w, value=E, gradient=dEdw)
}

