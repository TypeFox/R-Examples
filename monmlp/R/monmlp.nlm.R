monmlp.nlm <-
function(x, y, hidden1, hidden2=0, iter.max=5000, n.trials=1,
         Th=tansig, To=linear, Th.prime=tansig.prime, To.prime=linear.prime,
         monotone=NULL, init.weights=c(-0.5, 0.5), max.exceptions=10,
         silent=FALSE, ...)
{
    cost.best <- Inf
    for (i in 1:n.trials){
        exception <- TRUE
        n.exceptions <- 0
        while (exception){
            exception <- FALSE
            if (length(init.weights)==2){
                weights <- monmlp.initialize(x=x, y=y, hidden1=hidden1,
                                       hidden2=hidden2, init.weights=init.weights)
            } else{
                weights <- init.weights
            }
            output.mlp <- try(suppressWarnings(nlm(monmlp.cost, weights,
                              iterlim=iter.max, x=x, y=y, hidden1=hidden1,
                              hidden2=hidden2, Th=Th, To=To, Th.prime=Th.prime,
                              To.prime=To.prime, monotone=monotone,
                              check.analyticals=FALSE, ...)), silent=FALSE)
            if (class(output.mlp)=="try-error"){
                exception <- TRUE
                n.exceptions <- n.exceptions + 1
                if (n.exceptions > max.exceptions)
                    stop("maximum number of exceptions reached")
            }
        }
        weights <- output.mlp$estimate
        code <- output.mlp$code
        cost <- output.mlp$minimum
        if (!silent) cat(i, cost, "\n")
        if (cost < cost.best){
            weights.best <- weights
            cost.best <- cost
            code.best <- code
        }
    }
    list(weights=weights.best, cost=cost.best, code=code.best)
}

