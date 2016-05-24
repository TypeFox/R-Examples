enparCensored.km <-
function (x, censored, censoring.side, correct.se, left.censored.min, 
    right.censored.max, ci, ci.type, conf.level, pivot.statistic, 
    ci.sample.size = ifelse(pivot.statistic == "t", n.not.cen, 
        N)) 
{
    ppoints.list <- ppointsCensored(x = x, censored = censored, 
        censoring.side = censoring.side, prob.method = "kaplan-meier")
    ord.stats <- ppoints.list$Order.Statistics
    F.x <- ppoints.list$Cumulative.Probabilities
    Censored <- ppoints.list$Censored
    N <- length(ord.stats)
    n.cen <- sum(Censored)
    n.not.cen <- N - n.cen
    ord.stats.no.cen <- ord.stats[!Censored]
    F.x.no.cen <- F.x[!Censored]
    ord.stats.cen <- ord.stats[Censored]
    y <- ord.stats.no.cen
    F.y <- F.x.no.cen
    r <- (1:N)[!Censored][-n.not.cen]
    ci.ub <- Inf
    if (censoring.side == "left") {
        dl <- min(ord.stats.cen)
        if (dl <= min(ord.stats.no.cen) && left.censored.min != 
            "Ignore") {
            index <- match(dl, ord.stats.cen)
            F.y <- c(F.x[Censored][index], F.y)
            if (is.numeric(left.censored.min)) {
                dl <- left.censored.min
            }
            else {
                if (left.censored.min == "DL/2") 
                  dl <- dl/2
            }
            y <- c(dl, y)
            r <- c(1, r)
        }
    }
    else {
        dl <- max(ord.stats.cen)
        if (dl >= max(ord.stats.no.cen) && right.censored.max != 
            "Ignore") {
            if (is.numeric(right.censored.max)) 
                dl <- right.censored.max
            y <- c(y, dl)
            F.y <- c(F.y, 1)
            r <- c(r, (1:N)[!Censored][n.not.cen])
            ci.type <- "lower"
            ci.ub <- dl
        }
    }
    y <- c(0, y)
    F.y <- c(0, F.y)
    S.y <- 1 - F.y
    F.y.diffs <- diff(F.y)
    S.y.times.diff.y <- S.y[-length(S.y)] * diff(y)
    mu.hat <- sum(S.y.times.diff.y)
    sigma.hat <- sqrt(sum((y[-1] - mu.hat)^2 * F.y.diffs))
    A <- rev(cumsum(rev(S.y.times.diff.y[-1])))
    se.muhat <- sqrt(sum(A^2/((N - r) * (N - r + 1))))
    if (correct.se) 
        se.muhat <- sqrt(n.not.cen/(n.not.cen - 1)) * se.muhat
    parameters <- c(mean = mu.hat, sd = sigma.hat, se.mean = se.muhat)
    if (ci) {
        ci.obj <- ci.normal.approx(theta.hat = parameters["mean"], 
            sd.theta.hat = parameters["se.mean"], n = ci.sample.size, 
            df = ci.sample.size - 1, ci.type = ci.type, alpha = 1 - 
                conf.level, lb = 0, ub = ci.ub, test.statistic = pivot.statistic)
        ci.obj$parameter <- "mean"
        return(list(parameters = parameters, ci.obj = ci.obj))
    }
    else return(list(parameters = parameters))
}
