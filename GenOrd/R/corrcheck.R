corrcheck<-function (marginal, support = list(), Spearman = FALSE) 
{
    if (!all(unlist(lapply(marginal, function(x) (sort(x)==x & min(x)>0 & max(x)<1))))) stop("Error in assigning marginal distributions!")
    k <- length(marginal)
    mcmax <- diag(1, k)
    mcmin <- diag(1, k)
    len <- length(support)
    if (len == 0) {
        for (i in 1:k) {
            support[[i]] <- 1:(length(marginal[[i]]) + 1)
        }
    }
    if (Spearman) {
        for (i in 1:k) {
        	s1<-   c(marginal[[i]], 1)
		s2<-   c(0,marginal[[i]])
     		support[[i]] <- (s1+s2)/2
        }
    }
    for (i in 1:(k - 1)) {
        for (j in (i + 1):k) {
            P1 <- c(0, marginal[[i]], 1)
            P2 <- c(0, marginal[[j]], 1)
            l1 <- length(P1) - 1
            l2 <- length(P2) - 1
            p1 <- numeric(0)
            p2 <- numeric(0)
            for (g in 1:l1) {
                p1[g] <- P1[g + 1] - P1[g]
            }
            for (g in 1:l2) {
                p2[g] <- P2[g + 1] - P2[g]
            }
            E1 <- sum(p1 * support[[i]])
            E2 <- sum(p2 * support[[j]])
            V1 <- sum(p1 * support[[i]]^2) - E1^2
            V2 <- sum(p2 * support[[j]]^2) - E2^2
            y1 <- 1
            y2 <- 1
            lim <- 0
            E12 <- 0
            PP1 <- P1
            PP2 <- P2
            PP1 <- PP1[-1]
            PP2 <- PP2[-1]
            while (length(PP1) > 0) {
                E12 <- E12 + support[[i]][y1] * support[[j]][y2] * 
                  (min(PP1[1], PP2[1]) - lim)
                lim <- min(PP1, PP2)
                if (PP1[1] == lim) {
                  PP1 <- PP1[-1]
                  y1 <- y1 + 1
                }
                if (PP2[1] == lim) {
                  PP2 <- PP2[-1]
                  y2 <- y2 + 1
                }
            }
            c12 <- (E12 - E1 * E2)/sqrt(V1 * V2)
            y1 <- 1
            y2 <- l2
            lim <- 0
            E21 <- 0
            PP1 <- P1
            PP2 <- cumsum(rev(p2))
            PP1 <- PP1[-1]
            while (length(PP1) > 0) {
                E21 <- E21 + support[[i]][y1] * support[[j]][y2] * 
                  (min(PP1[1], PP2[1]) - lim)
                lim <- min(PP1, PP2)
                if (PP1[1] == lim) {
                  PP1 <- PP1[-1]
                  y1 <- y1 + 1
                }
                if (PP2[1] == lim) {
                  PP2 <- PP2[-1]
                  y2 <- y2 - 1
                }
            }
            c21 <- (E21 - E1 * E2)/sqrt(V1 * V2)
            mcmax[i, j] <- c12
            mcmin[i, j] <- c21
        }
    }
    mcmax <- forceSymmetric(mcmax)
    mcmin <- forceSymmetric(mcmin)
    list(mcmin, mcmax)
}
