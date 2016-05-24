`fblr` <-
function (y, bin, niter, thin = 5, nburn = 10000, int.level = 2, 
    kmax = 10, geo = 1, delta1 = 0.001, delta2 = 0.1, predict = FALSE, 
    file = "fblr_mcmc.txt") 
{
    require(MASS)     # library MASS for mvrnorm
    In <- function(n) diag(rep(1, n))
    # Generate z values from truncated normal with the Inverse cdf method
    trunc.norm <- function(mu) {
        mu + qnorm((pnorm(-mu) + runif(length(mu)) * (1 - pnorm(-mu))))
    }
    generate.z <- function(y, eta) {
        ifelse(y == 1, trunc.norm(eta), -trunc.norm(-eta))
    }
    
    
    # the non-Bayesfactor part of the acceptance probability
    compute.R <- function(model, mi, geo){
        R <- 1
        sizes <- model[mi$begin]
        if (model[1] == 0 | mi$model[1] == 0) { # Move to or from null model
            R <- switch(2 + 1 * (model[1] == 0) - 1 * (mi$model[1] == 
                0), 4, 1, 1/4)
        }
        else if (model[1] == mi$model[1]){      # change
            if (sum(sizes) < sum(mi$model[mi$begin])) 
                R <- geo                        # cbirth
            else if (sum(sizes) > sum(mi$model[mi$begin])) 
                R <- 1/geo                      # cdeath
        }
        else if (model[1] < mi$model[1]) 
            R <- 1/(sum(sizes == 1) + 1)        # birth
        else if (model[1] > mi$model[1]) 
            R <- sum(sizes==1)                  # death
        R
    }

    #Compute acceptance probability
    acceptance4 <- function(model.old, move.info, marg.loglike, 
        marg.loglike.new, geo) {
        min(c(1, compute.R(model.old, move.info, geo = geo) * 
            exp(marg.loglike.new - marg.loglike)))   # Bayes Factor
    }
    
    compare <- function(x, y) (all(x == y)) * 1
    eval.logic.single <- function(logic, bin) {
        size <- logic[1]
        if (size == 0) 
            return(NULL)
        lvars <- logic[2:(size + 1)]
        norm <- (lvars > 0) * 1
        if (size == 1) 
            return((bin[, abs(lvars)] == norm) * 1)
        else apply(bin[, abs(lvars)], 1, compare, y = norm)
    }
   
    # update model matrix
    change.matrix <- function(mi, B) {
        if (mi$model[1] == 0) 
            return(matrix(rep(1, n), ncol = 1))
        if (mi$move.type == 3) 
            B[, mi$pick + 1] <- eval.logic.single(mi$model[mi$begin[mi$pick]:(mi$begin[mi$pick] + int.level)], bin = bin)
        if(mi$move.type == 2) 
            B <- cbind(B, if (mi$k == 0) 
                eval.logic.single(mi$model[2:(2 + int.level)], bin = bin)
            else eval.logic.single(mi$model[(mi$begin[mi$k] + int.level + 
                1):(mi$begin[mi$k] + 2 * int.level + 1)], bin = bin))
        if (mi$move.type == 1) {
            B[, mi$pick + 1] <- B[, mi$k + 1]
            B <- B[, 1:mi$k]
        }
        B
    }

    # generate a candidate model
    gen.cand <- function(model){ 
        k <- model[1]
        if(k == 0){
                model[1:3] <- c(1, 1, sample(c(-1, 1), 1) * sample(vars, 1))
                return(list(model=model,k=0, begin=NULL,  move.type=2, pick=NULL))    
                }
        begin <- seq(from = 2, length.out = k, by = int.level + 1)
        sizes <- model[begin]
        death.poss <- sum(sizes == 1)
            u <- runif(1)
            if (u < (death.poss > 0) * 0.25) 
                move.type <- 1
            else if (u < (death.poss > 0) * 0.25 + (k < kmax) * 
                0.25) 
                move.type <- 2
            else move.type <- 3
      switch(move.type, death(model,k,begin,sizes,death.poss), birth(model,k,begin), 
                change(model, k, begin, sizes))  
    }
    
    # functions for the different moves
    death <- function(model,k, begin, sizes, death.poss){
        dead <- ifelse(death.poss == 1, which(sizes == 1), sample(which(sizes == 1), 1))
        model[begin[dead]:(begin[dead] + int.level)] <- model[begin[k]:(begin[k] + int.level)]
        model[begin[k]:(begin[k] + int.level)] <- rep(0, int.level + 1)
        model[1] <- model[1] - 1
        list(model=model, k=k,begin=begin, move.type=1, pick= dead)
    }
   
   birth <- function(model,k,begin){
        model[(begin[k] + int.level + 1):(begin[k] + int.level + 
            2)] <- c(1, sample(c(-1, 1), 1) * sample(vars, 1))
        model[1] <- model[1] + 1
        list(model=model,k=k,begin=begin, move.type=2, pick=NULL)
    }
   
    change <- function(model,k,begin,sizes){
        pick <- sample(1:k, 1)
        size <- sizes[pick]
        u <- runif(1)
        if (u < 0.3333 * (size > 1)) 
            cmove <- 1
        else if (u < 0.3333 * (size > 1) + (size != int.level) * 
            0.3333) 
            cmove <- 2
        else cmove <- 3
        logic <- model[begin[pick]:(begin[pick] + int.level)]
        model[begin[pick]:(begin[pick] + int.level)] <- switch(cmove, 
            cdeath(logic), cbirth(logic), cchange(logic))
        list(model=model, k=k,begin=begin,move.type=3,pick=pick)
    }
    sort.abs <- function(vek) vek[order(abs(vek))]
    cdeath <- function(logic) {
        size <- logic[1]
        pickvar <- sample(1:size, 1)
        logic[c(1, pickvar + 1, size + 1)] <- c(size - 1, logic[size + 
            1], 0)
        c(logic[1], sort.abs(logic[2:size]), rep(0, int.level - 
            size + 1))
    }
    cbirth <- function(logic) {
        size <- logic[1]
        lvars <- abs(logic[2:(size + 1)])
        c(size + 1, sort.abs(c(logic[2:(size + 1)], sample(c(-1, 
            1), 1) * sample(vars[-lvars], 1))), rep(0, int.level - 
            size - 1))
    }
    cchange <- function(logic) {
        size <- logic[1]
        cpick <- sample(1:size, 1)
        lvars <- abs(logic[(2:(size + 1))[-cpick]])
        logic[cpick + 1] <- ifelse(size == 1, sample(c(-1, 1), 
            1) * sample(vars, 1), sample(c(-1, 1), 1) * sample(vars[-lvars], 
            1))
        c(size, sort.abs(logic[2:(size + 1)]), rep(0, int.level - 
            size))
    }
    
    # delete files from previous runs
    write.table(file = file, x = NULL, col.names = FALSE) 
     
    # Initialize
    nbin <- dim(bin)[2]
    vars <- 1:dim(bin)[2]
    n <- length(y)
    if (predict) 
        pred <- rep(0, n)
    model <- rep(0, 1 + kmax * (int.level + 1))
    k <- 0
    B <- as.matrix(rep(1, n))
    b <- 0
    accept <- 0
    # MCMC run
    for (i in 1:(niter + nburn)) {
        eta <- B %*% b
        z <- generate.z(y, eta)
        tau <- rgamma(1, shape = delta1 + 0.5 * (k + 1), rate = delta2 + 
            0.5 * t(b) %*% b)
        Vstar <- solve(t(B) %*% B + tau * In(k + 1))
        bstar <- t(z) %*% z - t(z) %*% B %*% Vstar %*% t(B) %*% z
        marg.loglike <- 0.5 * (log(abs(det(Vstar))) - bstar + (k + 1) * log(tau))
        
        move.info <- gen.cand(model)
        B.new <- change.matrix(move.info, B = B)
        k.new <- dim(B.new)[2] - 1
        Vstar.new <- solve(t(B.new) %*% B.new + tau * In(k.new + 1))
        bstar.new <- t(z) %*% z - t(z) %*% B.new %*% Vstar.new %*% 
            t(B.new) %*% z
        marg.loglike.new <- 0.5 * (log(abs(det(Vstar.new))) - 
            bstar.new + (k.new + 1) * log(tau))
        if (runif(1) < acceptance4(model, move.info, marg.loglike, 
            marg.loglike.new, geo = geo)) {
            model <- move.info$model
            B <- B.new
            k <- k.new
            Vstar <- Vstar.new
            marg.loglike <- marg.loglike.new
            accept <- accept + 1
        }
        m <- Vstar %*% t(B) %*% z
        b <- mvrnorm(1, mu = m, Sigma = Vstar)
        if (i%%thin == 0 & i > nburn) {
            write.table(matrix(c(model, round(marg.loglike, 2), 
                round(tau, 6), round(b, 4), rep(0, kmax + 1 - 
                  length(b))), nrow = 1), file = file, append = TRUE, 
                row.names = FALSE, col.names = FALSE)
            if (predict) 
                pred <- pred + pnorm(B %*% b)
        }
    }
    if(predict)
	return(list(accept = accept/(niter + nburn), pred = as.vector(pred)/(niter/thin)))
    else
        return(accept/(niter + nburn))
}

