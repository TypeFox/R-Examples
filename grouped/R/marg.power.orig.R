"marg.power.orig" <-
function(MC, n, m, theta, sigma, a, dist.t, dist.x, grouping.mech){
    treat <- if(dist.t[1, 1] == "bernoulli") rbinom(n, prob = dist.t[1, 2], size = dist.t[1, 3]) else
             if(dist.t[1, 1] == "no distr") sample(c(rep(0, dist.t[1, 2]), rep(1, dist.t[1, 3])))
    if(missing(dist.x)) X <- cbind(rep(1, n), treat) else{
        n.covs <- nrow(dist.x)
        X <- matrix(0, ncol = n.covs, nrow = n)
        for(k in 1:n.covs){
            rdist <- if(dist.x[k, 1] == "normal")  "rnorm" else
            if(dist.x[k, 1] == "gamma") "rgamma" else
            if(dist.x[k, 1] == "beta") "rbeta" else
            if(dist.x[k, 1] == "chisquare") "rchisq" else 
            if(dist.x[k, 1] == "uniform") "runif" else
            if(dist.x[k, 1] == "bernoulli") "rbinom"
            X[, k] <- eval(parse(text = paste(rdist, "(n,", dist.x[k, 2], ",", dist.x[k, 3], ")", sep = "")))
        }
        X <- cbind(rep(1, n), treat, X)
    }
    power.. <- numeric(MC)  
    for(j in 1:MC){  
         power..[j] <- cond.power.orig(X = X, theta = theta, sigma = sigma, m = m, a = a, grouping.mech = grouping.mech)
    }
    list(mean(power..), power..)
}

