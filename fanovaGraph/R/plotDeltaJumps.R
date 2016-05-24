plotDeltaJumps <- function(graphlist, interval = c(0,1), mean.clique.size = FALSE) {
    op <- par(no.readonly = TRUE)
    par(mfrow = c(2, 1), mar = c(0, 4, 1, 0), oma = c(6, 0.5, 2, 1), 
        mgp = c(2.5, 1, 0))
    d <- graphlist$d
    dall <- graphlist$V
    totalInt <- rbind(combn(d,2),t(graphlist$tii[,1]))
    o <- order(totalInt[3, ], decreasing = FALSE)
    totalIntSort <- totalInt[, o]
    totalIntSortNorm <- rbind(totalIntSort[1:2, ], totalIntSort[3, 
        ]/dall)
    delta <- c(0, totalIntSortNorm[3, ], 1)
    delta <- delta[delta >= interval[1] & delta <= interval[2]]
    n.CL <- c()
    s.CL <- c()
    for (i in 1:length(delta)) {
        E <- t(totalIntSortNorm[-3, which(totalIntSortNorm[3, ] > delta[i])])
        E.graph <- graph(as.vector(t(E)), n = d , directed = FALSE)
        CL <- maximal.cliques(E.graph)
        n.CL[i] <- length(CL)
        s.CL[i] <- mean(sapply(CL, length))
    }
    plot(1:length(delta), n.CL, type = "s", ylab = "number of cliques", 
        xaxt = "n", ylim = c(0, max(n.CL)))
    if (mean.clique.size == TRUE) {
        lines(1:length(delta), s.CL, type = "s", lty = 3, lwd = 2)
        legend("topright", c("number of cliques", "mean clique size"), 
            lty = c(1, 3))
    }
    deltaCut <- c(delta[-length(delta)], delta[length(delta) - 1])
    jumps <- diff(deltaCut)
    alphas <- jumps/max(jumps)*0.3
    alphas <- sqrt(alphas)    
    for (i in 1:(length(delta)-1))
      polygon(c(i,i+1,i+1,i), c(-1,-1,max(n.CL)*10,max(n.CL)*10), 
              col = rgb(0,1,0,alpha=alphas[i]), border = NA)
    
    plot(1:length(delta), deltaCut, type = "s", ylab = "jumps", xaxt = "n")
    axis(1, at = 1:length(delta), labels = round(delta, 4), las = 2, 
        outer = TRUE)
    title(xlab = "steps in delta", outer = TRUE, mgp = c(4, 1, 0))
    title("Delta Jump Plot", outer = TRUE)
    for (i in 1:(length(delta)-1))
      polygon(c(i,i+1,i+1,i), c(-1,-1,max(n.CL)*10,max(n.CL)*10), 
              col = rgb(0,1,0,alpha=alphas[i]), border = NA)
    par(op)
} 