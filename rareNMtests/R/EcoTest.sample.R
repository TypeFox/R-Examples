EcoTest.sample <-
function (x, by = NULL, MARGIN = 2, niter = 200, method = "sample-size", q = 0, trace = TRUE) {

    if (is.null(by))
        stop("grouping of samples is not specified in 'by'")
    if (length(by) != dim(x)[MARGIN])
        stop("length of 'by' does not match number of sites in matrix")

    if (MARGIN == 1) {
        x <- as.matrix(x[order(by),])
    } else {
        x <- as.matrix(t(x[,order(by)]))
    }
    x <- ifelse(x>0, 1, 0)

    obs <- list()
    for(i in 1:length(levels(factor(by)))) {
        obs[[i]] <- rarefaction.sample(subset(x, by==levels(factor(by))[i]), method = method, q = q)
        }
    if (method == "coverage") {
        obs <- lapply(obs, function(x) rbind(c(0,0), x))
    }
    obs.t <- rarefaction.sample(x, method = method, q = q)
    if (method == "coverage") {
        obs.t <- rbind(c(0,0), obs.t)
    }

    A <- NULL
    for(i in 1:length(levels(factor(by)))) {
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

    Zsim <- NULL
    for (i in 1:niter) {
        if (trace == TRUE) 
            print(paste("This is randomization number", i))
        random.m <- x[sample(1:dim(x)[1], dim(x)[1]), ]
        sim <- list()
        for(z in 1:length(levels(factor(by)))) {
            sim[[z]] <- rarefaction.sample(subset(random.m, by==levels(factor(by))[z]), method = method, q = q)
            }
        if (method == "coverage") {
            sim <- lapply(sim, function(x) rbind(c(0,0), x))
        }

        A <- NULL
        for(j in 1:length(levels(factor(by)))) {
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
        M <- "Sample-based method"
        cat(M, "\n")
        cat("P(Obs <= null) < ", format(p.upper, digits=5), "\n")
    } else if (Z < min(Zsim)) {
        p.lower <- 1/length(Zsim) 
        p.upper <- (length(Zsim) - 1)/length(Zsim)
        M <- "Sample-based method"
        cat(M, "\n")
        cat("P(Obs <= null) > ", format(p.upper, digits=5), "\n")
    } else {
        p.lower <- sum(Z >= Zsim)/length(Zsim)
        p.upper <- sum(Z <= Zsim)/length(Zsim)
        M <- "Sample-based method"
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
