BiogTest.sample <-
function (x, by = NULL, MARGIN = 2, niter = 200, method = "sample-size", q = 0, trace = TRUE, distr = "lnorm") {

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
    A <- NULL
    p <- combn(1:length(obs), m = 2)
    lrc <- do.call("c", lapply(obs, function(t) max(t[,1])))
    for(i in 1:dim(p)[2]) {
        max.cov <- min(lrc[p[,i]])
        nxax1 <- findInterval(max.cov, obs[[p[1,i]]][,1])
        nxax2 <- findInterval(max.cov, obs[[p[2,i]]][,1])
        xax1 <- obs[[p[1,i]]][1:nxax1,1]
        xax2 <- obs[[p[2,i]]][1:nxax2,1]
        yval1 <- obs[[p[1,i]]][1:nxax1,2]
        yval2 <- obs[[p[2,i]]][1:nxax2,2]
        poly <- cbind(c(xax2, rev(xax1)), c(yval2, rev(yval1)))
        A <- c(A, abs(area.xypolygon(poly)))
        }
    Z <- sum(A)

    SR <- NULL
    for(i in 1:length(levels(factor(by)))) {
        SR <- rbind(SR, unlist(specpool(subset(x, by==levels(factor(by))[i]))[2:3]))
    }
    mr <- match(max(SR[,1]), SR[,1])
    range.SR <- c(SR[mr, 1]-2*SR[mr,2], SR[mr, 1]+2*SR[mr,2])
    obs.SR <- by(x, factor(as.character(by)), function(z) ncol(z[,colSums(z)>0]))
    range.SR[1] <- ifelse(range.SR[1] < max(obs.SR), max(obs.SR), range.SR[1]) 

    Zsim <- NULL
    for (it in 1:niter) {
        if (trace == TRUE) 
            print(paste("This is randomization number", it))
    repeat {
        richness <- round(runif(1, range.SR[1], range.SR[2]))
        if (distr == "lnorm") {
            sdlog <- runif(1, 0.1, 3.5)
            com <- rlnorm(richness, 0, sdlog)
        } else
        if (distr == "geom") {
            p <- runif(1, 0.1, 1)
            com = NULL
            com[1]= p
            for (i in 2:richness)
                com[i] = com[i - 1] * p
        } else 
        if (distr == "stick") {
            segments <- runif(richness-1)
            segments <- c(0, segments, 1)
            segments <- sort(segments)
            com <- diff(segments)
        }
        com <- ceiling(com/min(com[com>0]))
        prob <- rnbinom(n = length(com), mu=com, size=runif(length(com), 0.1, 25))
        index <- is.na(prob)
        df <- sample(1:length(com[!index]), sum(x[1,]), replace=FALSE, prob=prob[!index])
        df <- data.frame(table(df))
        colnames(df) <- c("species", "sample1")
        for(j in 2:dim(x)[1]) {
            prob <- rnbinom(n = length(com), mu=com, size=runif(length(com), 0.1, 25))
            index <- is.na(prob)
            stemp <- sample(1:length(com[!index]), sum(x[j,]), replace=FALSE, prob=prob[!index])
            stemp <- data.frame(table(stemp))
            colnames(stemp) <- c("species", paste("sample", j, sep=""))
            df <- merge(df, stemp, by="species", all=TRUE)
        }
        df[is.na(df)] <- 0
        df <- t(df[,-1])
        dfv <- aggregate(df>0, list(by), sum)
        dfv <- apply(dfv[,-1]>0, 1, sum)
        if(any(dfv == 1)) {
            next
            }
        sim <- list()
        for(i in 1:length(levels(factor(by)))) {
            sim[[i]] <- rarefaction.sample(subset(df, by==levels(factor(by))[i]), method = method, q = q)
            }
        simv <- sum(unlist(lapply(sim, function(x) sum(diff(order(x[,1]))!=1))))
        if(simv !=0) {
            next
            }
        A <- NULL
        p <- combn(1:length(sim), m = 2)
        lrc <- do.call("c", lapply(obs, function(t) max(t[,1])))
        for(k in 1:dim(p)[2]) {
            max.cov <- min(lrc[p[,k]])
            nxax1 <- findInterval(max.cov, sim[[p[1,k]]][,1])
            nxax2 <- findInterval(max.cov, sim[[p[2,k]]][,1])
            xax1 <- sim[[p[1,k]]][1:nxax1,1]
            xax2 <- sim[[p[2,k]]][1:nxax2,1]
            yval1 <- sim[[p[1,k]]][1:nxax1,2]
            yval2 <- sim[[p[2,k]]][1:nxax2,2]
            poly <- cbind(c(xax2, rev(xax1)), c(yval2, rev(yval1)))
            A <- c(A, abs(area.xypolygon(poly)))
            }
        Zsim <- c(Zsim, sum(A))
        if(length(Zsim)==it) {break}
    }
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

    output <- list(subclass = M, type = k, obs = obs, sim = sim, Zsim = Zsim, Z = Z, pval = p.values)
    class(output) <- c("BiogTest", class(output))
    return(output)
}
