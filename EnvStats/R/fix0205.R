fix0205 <-
function (data, censor, alpha, beta = 1) 
{
    start <- date()
    estimate <- matrix(-1, 1, 6)
    mestimate <- matrix(-1, 1, 6)
    realvalues <- data[data > censor]
    pqlvalues <- data[data <= censor]
    if (length(realvalues) > 0 & length(pqlvalues) > 0) 
        tester <- 1
    if (length(realvalues) == 0 & length(pqlvalues) > 0) 
        tester <- 2
    if (length(realvalues) > 0 & length(pqlvalues) == 0) 
        tester <- 3
    allestim <- switch(tester, fix0101(cbind(data, rep(censor, 
        length(data))), alpha = alpha, beta = beta), stop("All values are censored"), 
        fix0201.obs(cbind(data, rep(censor, length(data))), alpha = alpha, 
            beta = beta))
    estimate[1, 1:2] <- c(alpha, beta)
    estimate[1, 3:6] <- allestim$conestimate
    mestimate[1, 1:6] <- allestim$mestimate
    labl <- c("alpha", "beta", "mean", "std.dev.", "exp.psi", 
        "exp.chi")
    dimnames(estimate) <- list(1, labl)
    mlabl <- c("obs", "iter", "mu", "sigma", "psi", "chi")
    dimnames(mestimate) <- list(1, mlabl)
    sample <- c(n = sum(!is.na(data)), nmiss = sum(is.na(data)), 
        censored = sum(data <= censor), uncensored = sum(data > 
            censor), censor.val = censor)
    summary <- summary(data)
    asympt.var <- fix0202(mu = mestimate[1, 3], sig = mestimate[1, 
        4], eta = estimate[1, 3], kappa = estimate[1, 4], alf = alpha, 
        bet = beta, mdl = censor)
    consistent.deriv.i <- solve(fix0207(mu = mestimate[1, 3], 
        sig = mestimate[1, 4], eta = estimate[1, 3], kappa = estimate[1, 
            4], alf = alpha, bet = beta, mdl = censor))
    asympt.consistent <- consistent.deriv.i %*% asympt.var$psi.var %*% 
        t(consistent.deriv.i)
    dimnames(asympt.var$asympt.var) <- list(c("mu", "sig"), c("mu", 
        "sig"))
    dimnames(asympt.consistent) <- list(c("mean", "std.dev."), 
        c("mean", "std.dev."))
    stderr.mu <- sqrt(asympt.var$asympt.var[1, 1]/sample[1])
    stderr.sig <- sqrt(asympt.var$asympt.var[2, 2]/sample[1])
    stderr.eta <- sqrt(asympt.consistent[1, 1]/sample[1])
    stderr.kappa <- sqrt(asympt.consistent[2, 2]/sample[1])
    list(t.df = 2 * alpha, sample = sample, consistent.estimates = estimate[, 
        3:4], covariance.estimates = asympt.consistent/length(data), 
        convergence.diag = estimate[, 5:6])
}
