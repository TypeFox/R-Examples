EcoTest.individual <-
function (x, MARGIN = 2, niter = 200, method = "sample-size", q = 0, trace = TRUE, powerfun = 1, log.scale = FALSE) {

    if ((sum(x>0) - sum(x)) == 0) 
        stop("matrix must contain abundances not presence-absence")
    if (MARGIN == 1) {
        x <- as.matrix(t(x))
    } else {
        x <- as.matrix(x)
    }

    obs <- apply(x, 2, function(x) rarefaction.individual(x, method = method, q = q, powerfun = powerfun, log.scale = log.scale))
    if (method == "coverage") {
        obs <- lapply(obs, function(x) rbind(c(0,0), x))
    }
    obs.t <- rarefaction.individual(rowSums(x), method = method, q = q, powerfun = powerfun, log.scale = log.scale)
    if (method == "coverage") {
        obs.t <- rbind(c(0,0), obs.t)
    }

    A <- NULL
    for(i in 1:ncol(x)) {
        max.cov <- min(max(obs[[i]][,1]), max(obs.t[,1]))
        nxax <- findInterval(max.cov, obs[[i]][,1])
        nxaxt <- findInterval(max.cov, obs.t[,1])
        xax <- obs[[i]][1:nxax,1]
        xaxt <- obs.t[1:nxaxt,1]
        yval1 <- obs[[i]][1:nxax,2]
        yval2 <- obs.t[1:nxaxt,2]
        poly <- cbind(c(xaxt, rev(xax)), c(yval2, rev(yval1)))
        A <- c(A, abs(area.xypolygon(poly)))
    }
    Z <- sum(A)

    indsSite <- colSums(x)
    indsSp <- rowSums(x)
    sites <- rep(indsSite, times=indsSite)
    sp <- rep(indsSp, times=indsSp)
    Zsim <- NULL
    for (i in 1:niter) {
        if (trace == TRUE) 
            print(paste("This is randomization number", i))
        random.m <- xtabs(~sample(names(sp), length(sp)) + sites)
        sim <- apply(random.m, 2, function(rm) rarefaction.individual(rm, method = method, q = q, powerfun = powerfun, log.scale = log.scale))
        if (method == "coverage") {
            sim <- lapply(sim, function(x) rbind(c(0,0), x))
        }

        A <- NULL
        for(j in 1:ncol(random.m)) {
            max.cov <- min(max(sim[[j]][,1]), max(obs.t[,1]))
            nxax <- findInterval(max.cov, sim[[j]][,1])
            nxaxt <- findInterval(max.cov, obs.t[,1])
            xax <- sim[[j]][1:nxax,1]
            xaxt <- obs.t[1:nxaxt,1]
            yval1 <- sim[[j]][1:nxax,2]
            yval2 <- obs.t[1:nxaxt,2]
            poly <- cbind(c(xaxt, rev(xax)), c(yval2, rev(yval1)))
            A <- c(A, abs(area.xypolygon(poly)))
        }
    Zsim <- c(Zsim, sum(A))
    }

    if (Z > max(Zsim)) {
        p.lower <- (length(Zsim) - 1)/length(Zsim) 
        p.upper <- 1/length(Zsim)
        M <- "Individual-based method"
        cat(M, "\n")
        cat("P(Obs <= null) < ", format(p.upper, digits=5), "\n")
    } else if (Z < min(Zsim)) {
        p.lower <- 1/length(Zsim) 
        p.upper <- (length(Zsim) - 1)/length(Zsim)
        M <- "Individual-based method"
        cat(M, "\n")
        cat("P(Obs <= null) > ", format(p.upper, digits=5), "\n")
    } else {
        p.lower <- sum(Z >= Zsim)/length(Zsim)
        p.upper <- sum(Z <= Zsim)/length(Zsim)
        M <- "Individual-based method"
        cat(M, "\n")
        cat("P(Obs <= null) = ", format(p.upper, digits=5), "\n")
    }
    p.values <- c(p.lower, p.upper)
    names(p.values) <- c("lower p_val", "upper p_val")
    k <- ifelse(method == "coverage", "Coverage measure", "Sample-size")

    output <- list(subclass = M, type = k, obs = obs, sim = sim, pooled = obs.t, Zsim = Zsim, Z = Z, pval = p.values)
    class(output) <- c("EcoTest", class(output))
    return(output)
}
