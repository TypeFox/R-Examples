gevcdn.bootstrap <-
function (n.bootstrap, x, y, iter.max = 1000, n.hidden = 2,
          Th = gevcdn.logistic, fixed = NULL,
          init.range = c(-0.25, 0.25), scale.min = .Machine$double.eps,
          beta.p = 3.3, beta.q = 2, sd.norm = Inf, n.trials = 5,
          method = c("BFGS", "Nelder-Mead"),
          boot.method = c("residual", "parametric"),
          init.from.prev = TRUE, max.fails = 100,
          probs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
          ...)
{
    if (!is.matrix(x)) stop("\"x\" must be a matrix")
    if (!is.matrix(y)) stop("\"y\" must be a matrix")
    boot.method <- match.arg(boot.method)
    weights.bootstrap <- parms.bootstrap <- quantiles.bootstrap <- list()
    location.bootstrap <- scale.bootstrap <- shape.bootstrap <- c()
    for (i in seq_len(n.bootstrap)){
        cat("** Trial", i, "\n")
        if (i==1){
            cat("Fitting model...\n")
            weights <- gevcdn.fit(x, y, iter.max, n.hidden, Th, fixed,
                                  init.range, scale.min,
                                  beta.p, beta.q, sd.norm,
                                  n.trials, method, max.fails, ...)
            parms <- gevcdn.evaluate(x, weights)
            residuals <- (1 + parms[,"shape"]*(y - parms[,"location"])/
                          parms[,"scale"])^(-1/parms[,"shape"])
        }
        if (boot.method=="residual"){
            y.prime <- as.matrix(parms[,"location"] + parms[,"scale"]*
                                (sample(residuals, replace=TRUE)^
                                (-parms[,"shape"]) - 1)/parms[,"shape"])
        } else if (boot.method=="parametric"){
            y.prime <- y*0
            for(j in seq_len(nrow(y))){
                  y.prime[j] <- rgev(1, location = parms[j,"location"],
                                        scale = parms[j,"scale"],
                                        shape = parms[j,"shape"])
            }
        }
        if (init.from.prev){
            n.trials <- 1
            init.range <- weights
        }
        weights.prime <- gevcdn.fit(x, y.prime, iter.max, n.hidden, Th,
                                    fixed, init.range, scale.min,
                                    beta.p, beta.q, sd.norm,
                                    n.trials, method, max.fails, ...)
        parms.prime <- gevcdn.evaluate(x, weights.prime)
        quantiles.prime <- sapply(probs, qgev,
                                  location=parms.prime[,"location"],
                                  scale=parms.prime[,"scale"],
                                  shape=parms.prime[,"shape"])
        colnames(quantiles.prime) <- probs
        quantiles.bootstrap[[i]] <- quantiles.prime
        weights.bootstrap[[i]] <- weights.prime
        parms.bootstrap[[i]] <- parms.prime
        location.bootstrap <- cbind(location.bootstrap,
                                    parms.prime[,"location"])
        scale.bootstrap <- cbind(scale.bootstrap, parms.prime[,"scale"])
        shape.bootstrap <- cbind(shape.bootstrap, parms.prime[,"shape"])
    }
    list(weights.bootstrap = weights.bootstrap,
         parms.bootstrap = parms.bootstrap,
         location.bootstrap = location.bootstrap,
         scale.bootstrap = scale.bootstrap,
         shape.bootstrap = shape.bootstrap,
         quantiles.bootstrap = quantiles.bootstrap)
}

