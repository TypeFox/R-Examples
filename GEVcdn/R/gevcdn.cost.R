gevcdn.cost <-
function (weights, x, y, n.hidden, Th, fixed, scale.min, beta.p, beta.q,
          sd.norm)
{
    w <- gevcdn.reshape(x, weights, n.hidden)
    attr(w, "Th") <- Th
    attr(w, "fixed") <- fixed
    attr(w, "scale.min") <- scale.min
    attr(w, "x.min") <- rep(0, ncol(x))
    attr(w, "x.max") <- rep(1, ncol(x))
    attr(w, "y.min") <- 0
    attr(w, "y.max") <- 1
    params <- gevcdn.evaluate(x, w)
    location <- params[,"location"]
    scale <- params[,"scale"]
    shape <- params[,"shape"]
    L <- dgev(y, location = location, scale = scale, shape = shape)
    if (!is.null(c(beta.p, beta.q))){
        prior <- dbeta(shape + 0.5, shape1 = beta.p, shape2 = beta.q)
        penalty <- -mean(log(prior))

    } else{
        penalty <- 0
    }
    if (sd.norm != Inf){
        prior.W1 <- dnorm(as.vector(w$W1), sd = sd.norm)
        penalty.W1 <- -mean(log(prior.W1))
    } else{
        penalty.W1 <- 0
    }
    NLL <- -sum(log(L))
    if(is.nan(NLL)) NLL <- .Machine$double.xmax
    GML <- NLL + penalty + penalty.W1
    attr(GML, "NLL") <- NLL
    attr(GML, "penalty") <- penalty + penalty.W1
    GML
}

