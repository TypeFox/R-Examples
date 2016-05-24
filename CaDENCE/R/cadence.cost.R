cadence.cost <-
function(weights, x, y, n.hidden, hidden.fcn, distribution, sd.norm, valid)
{
    weights.valid <- valid*0.
    weights.valid[valid] <- weights
    w <- cadence.reshape(x, weights.valid, n.hidden, distribution)
    cdn <- cadence.evaluate(x, w$W1, w$W2, hidden.fcn, distribution)
    args <- as.list(data.frame(cbind(y, cdn)))
    names(args) <- NULL
    L <- do.call(distribution$density.fcn, args)
    NLL <- -sum(log(L))
    penalty <- 0
    if(sd.norm != Inf){
        if(identical(hidden.fcn, identity)){
            penalty <- -mean(log(c(dnorm(as.vector(w$W1), sd=sd.norm),
                                   dnorm(as.vector(w$W2[1:(nrow(w$W2)-1),]),
                                         sd=sd.norm))))
        } else{
            penalty <- -mean(log(dnorm(as.vector(w$W1), sd=sd.norm)))
        }
    }
    if(is.nan(NLL)) NLL <- .Machine$double.xmax
    NLL <- NLL + penalty
    attr(NLL, "penalty") <- penalty
    NLL
}

