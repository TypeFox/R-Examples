`brier` <-
function (obs, pred, baseline = NULL, thresholds = seq(0, 1, 
    0.1), bins = TRUE, ...) 
{
    id <- is.finite(obs) & is.finite(pred)
    obs <- obs[id]
    pred <- pred[id]
    pred <- round(pred, 8)
    thresholds <- round(thresholds, 8 )
    if (max(pred) > 1 | min(pred) < 0) {
        cat("Predictions outside [0,1] range.  \n Are you certain this is a probability forecast? \n")
    }
    if (is.null(baseline)) {
        obar <- mean(obs)
        baseline.tf <- FALSE
    }
    else {
        obar <- baseline
        baseline.tf <- TRUE
    }
    bs.baseline <- mean((obar - obs)^2)
    if (bins) {
        XX <- probcont2disc(pred, bins = thresholds)
        pred <- XX$new
        new.mids <- XX$mids
    }
    else {
        if (length(unique(pred)) > 20) {
            warning("More than 20 unique probabilities. This could take awhile.")
        }
    }
    N.pred <- aggregate(pred, by = list(pred), length)
    N.obs <- aggregate(obs, by = list(pred), sum)
    if (bins) {
        XX <- data.frame(Group.1 = new.mids, zz = rep(0, length(thresholds) - 
            1))
        XX$Group.1 <- as.factor(XX$Group.1)
        N.pred$Group.1 <- as.factor(N.pred$Group.1)
        N.obs$Group.1 <- as.factor(N.obs$Group.1)
        N.pred <- merge(XX, N.pred, all.x = TRUE)
        N.obs <- merge(XX, N.obs, all.x = TRUE)
    }
    else {
        XX <- data.frame(Group.1 = thresholds, zz = rep(0, length(thresholds)))
        XX$Group.1 <- as.factor(XX$Group.1)
        N.pred$Group.1 <- as.factor(N.pred$Group.1)
        N.obs$Group.1 <- as.factor(N.obs$Group.1)
        N.pred <- merge(XX, N.pred, all.x = TRUE)
        N.obs <- merge(XX, N.obs, all.x = TRUE)
    }
    obar.i <- N.obs$x/N.pred$x
    y.i <- as.numeric(as.character(N.obs$Group.1))
    bs <- mean((pred - obs)^2)
    n <- length(obs)
    ss <- 1 - bs/bs.baseline
    bs.rel <- sum(N.pred$x * (y.i - obar.i)^2, na.rm = TRUE)/n
    bs.res <- sum(N.pred$x * (obar.i - obar)^2, na.rm = TRUE)/n
    bs.uncert <- obar * (1 - obar)
    check <- bs.rel - bs.res + bs.uncert
    prob.y <- N.pred$x/n
    return(list(baseline.tf = baseline.tf, bs = bs, bs.baseline = bs.baseline, 
        ss = ss, bs.reliability = bs.rel, bs.resol = bs.res, 
        bs.uncert = bs.uncert, y.i = y.i, obar.i = obar.i, prob.y = prob.y, 
        obar = obar, thres = thresholds, check = check, bins = bins))
}
