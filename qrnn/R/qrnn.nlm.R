qrnn.nlm <-
function(x, y, n.hidden, tau, iter.max, n.trials, bag, lower, eps.seq,
         Th, Th.prime, penalty, trace, ...)
{
    cases <- 1:nrow(x)
    if (bag) cases <- sample(nrow(x), replace=TRUE)
    x <- x[cases,,drop=FALSE]
    y <- y[cases,,drop=FALSE]
    if(length(lower) > 1) lower <- lower[cases]
	eps.seq <- sort(eps.seq, decreasing=TRUE)
    cost.best <- Inf
    for(i in 1:n.trials){
        weights <- qrnn.initialize(x, y, n.hidden)
        if(any(lower != -Inf)){
            for(eps in eps.seq){
                fit <- nlm(qrnn.cost, weights, iterlim=iter.max,
                           x=x, y=y, n.hidden=n.hidden, tau=tau,
                           lower=-Inf, eps=eps, Th=Th,
                           Th.prime=Th.prime, penalty=penalty,
                           check.analyticals=FALSE, ...)
                weights <- fit$estimate
            }
        } 
        for(eps in eps.seq){
            fit <- nlm(qrnn.cost, weights, iterlim=iter.max,
                       x=x, y=y, n.hidden=n.hidden, tau=tau,
                       lower=lower, eps=eps, Th=Th,
                       Th.prime=Th.prime, penalty=penalty,
                       check.analyticals=FALSE, ...)
            weights <- fit$estimate
        }
        cost <- fit$minimum
        if(trace) cat(i, cost, "\n")
        if(cost < cost.best){
            cost.best <- cost
            weights.best <- fit$estimate
        }
    }
    if(trace) cat("*", cost.best, "\n")
    weights.best <- qrnn.reshape(x, y, weights.best, n.hidden)
    weights.best
}

