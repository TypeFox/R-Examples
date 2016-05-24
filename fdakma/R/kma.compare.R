kma.compare <-
function (x, y0 = NULL, y1 = NULL, n.clust = c(1, 2), warping.method = c("NOalignment", 
    "shift", "dilation", "affine"), similarity.method = "d1.pearson", 
    center.method = "k-means", seeds = NULL, optim.method = "L-BFGS-B", 
    span = 0.15, t.max = 0.1, m.max = 0.1, n.out = NULL, tol = 0.01, 
    fence = TRUE, iter.max = 100, show.iter = 0, plot.graph = 0, 
    nstart = 2, return.all = FALSE) 
{
    if (length(n.clust) == 0) {
        stop("length of n.clust must be positive")
    }
    for (ind in 1:length(n.clust)) {
        if (n.clust[ind] <= 0) {
            stop("n.clust must be positive")
        }
        if (n.clust[ind] != floor(n.clust[ind]) | n.clust[ind] != 
            ceiling(n.clust[ind])) {
            warning("the cluster number must be integer, the value has been approximated to the nearest integer to continue")
            n.clust[ind] <- round(n.clust[ind])
        }
    }
    if (return.all == TRUE) {
        warning("return.all assigned to FALSE, only the best result in each case will be considered")
    }
    n.clust <- sort(n.clust)
    sim <- NULL
    sim_NOalignment <- NULL
    sim_ONLYSHI <- NULL
    sim_ONLYDIL <- NULL
    Result.NOalignment <- NULL
    Result.shift <- NULL
    Result.dilation <- NULL
    Result.affine <- NULL
    if (!"NOalignment" %in% warping.method && !"shift" %in% warping.method && 
        !"dilation" %in% warping.method && !"affine" %in% warping.method) {
        stop("warping.method does not contain any accepted values")
    }
    cycle <- NULL
    colors <- NULL
    if ("NOalignment" %in% warping.method) {
        cycle <- c(cycle, "NOalignment")
        colors <- c(colors, "black")
    }
    if ("shift" %in% warping.method) {
        cycle <- c(cycle, "shift")
        colors <- c(colors, "blue")
    }
    if ("dilation" %in% warping.method) {
        cycle <- c(cycle, "dilation")
        colors <- c(colors, "forestgreen")
    }
    if ("affine" %in% warping.method) {
        cycle <- c(cycle, "affine")
        colors <- c(colors, "orange")
    }
    if (length(dim(y0)) != 0) {
        n.obs <- dim(y0)[1]
    }
    if (length(dim(y1)) != 0) {
        n.obs <- dim(y1)[1]
    }
    if (is.null(seeds)) {
        seeds <- sample(1:n.obs, max(n.clust))
    }
    if (length(dim(seeds)) == 0) {
        seeds <- as.matrix(t(seeds))
    }
    if (length(seeds[1, ]) > max(n.clust)) {
        stop("number of columns of \"seeds\" must be inferior or equal to max(n.clust)")
    }
    if (dim(seeds)[1] > nstart) {
        warning("Number of row of seeds higher than nstart, only the first nstart rows of seeds will be considered")
        seeds <- seeds[1:nstart, ]
        if (length(dim(seeds)) == 0) {
            seeds <- as.matrix(t(seeds))
        }
    }
    if (length(seeds[1, ]) < max(n.clust)) {
        seeds_2 <- matrix(NA, nrow = nstart, ncol = max(n.clust))
        seeds_2[1:dim(seeds)[1], 1:dim(seeds)[2]] <- seeds
        seeds_2[which(is.na(seeds_2))] <- sample(1:n.obs, length(which(is.na(seeds_2))))
        seeds <- seeds_2
    }
    if (length(seeds) != 0) {
        for (contr in 1:length(seeds)) {
            if (seeds[contr] > n.obs || seeds[contr] <= 0) {
                stop("At least a value of \"seeds\" is not valid (is negative or null or superior to the number of observations)")
            }
        }
    }
    optim.method.available <- c("L-BFGS-B", "SANN")
    if (length(optim.method) != 0) {
        if (!optim.method %in% optim.method.available) {
            stop("Value of \"optim.method\" not valid. If defined, it must be one of the following methods (optim package methods): ", 
                "\"", optim.method.available[1], "\"", " ", "\"", 
                optim.method.available[2], "\"")
        }
    }
    if (length(optim.method) == 0) {
        optim.method <- "L-BFGS-B"
    }
    similarity.method.compare <- similarity.method
    for (k in n.clust) {
        for (warping.method.tra in cycle) {
            return.all <- FALSE
            Result <- kma(x = x, y0 = y0, y1 = y1, n.clust = k, 
                warping.method = warping.method.tra, similarity.method = similarity.method.compare, 
                center.method = center.method, seeds = seeds[, 
                  1:k], optim.method = optim.method, span = span, 
                t.max = t.max, m.max = m.max, n.out = n.out, 
                tol = tol, fence = fence, iter.max = iter.max, 
                show.iter = show.iter, nstart = nstart, return.all = return.all)
            if (warping.method.tra == "NOalignment") {
                sim_NOalignment <- c(sim_NOalignment, mean(Result$similarity.final))
                Result.NOalignment <- c(Result.NOalignment, list(Result))
            }
            if (warping.method.tra == "shift") {
                sim_ONLYSHI <- c(sim_ONLYSHI, mean(Result$similarity.final))
                Result.shift <- c(Result.shift, list(Result))
            }
            if (warping.method.tra == "dilation") {
                sim_ONLYDIL <- c(sim_ONLYDIL, mean(Result$similarity.final))
                Result.dilation <- c(Result.dilation, list(Result))
            }
            if (warping.method.tra == "affine") {
                sim <- c(sim, mean(Result$similarity.final))
                Result.affine <- c(Result.affine, list(Result))
            }
        }
    }
    if (plot.graph == 1 && length(n.clust) > 1) {
        dev.new()
        massimo <- max(sim_NOalignment, sim_ONLYSHI, sim_ONLYDIL, 
            sim)
        minimo <- min(sim_NOalignment, sim_ONLYSHI, sim_ONLYDIL, 
            sim)
        for (indice in 1:length(cycle)) {
            if (indice == 1) {
                if (cycle[indice] == "NOalignment") 
                  plot(sim_NOalignment, type = "b", axes = F, 
                    ylim = c(minimo, massimo), xlab = "", ylab = "", 
                    col = "black", lwd = 4)
                if (cycle[indice] == "shift") 
                  plot(sim_ONLYSHI, type = "b", axes = F, ylim = c(minimo, 
                    massimo), xlab = "", ylab = "", col = "blue", 
                    lwd = 4)
                if (cycle[indice] == "dilation") 
                  plot(sim_ONLYDIL, type = "b", axes = F, ylim = c(minimo, 
                    massimo), xlab = "", ylab = "", col = "forestgreen", 
                    lwd = 4)
                if (cycle[indice] == "affine") 
                  plot(sim, type = "b", axes = F, ylim = c(minimo, 
                    massimo), xlab = "", ylab = "", col = "orange", 
                    lwd = 4)
            }
            else {
                if (cycle[indice] == "NOalignment") 
                  points(sim_NOalignment, type = "b", ylim = c(minimo, 
                    massimo), xlab = "", ylab = "", col = "black", 
                    lwd = 4)
                if (cycle[indice] == "shift") 
                  points(sim_ONLYSHI, type = "b", ylim = c(minimo, 
                    massimo), xlab = "", ylab = "", col = "blue", 
                    lwd = 4)
                if (cycle[indice] == "dilation") 
                  points(sim_ONLYDIL, type = "b", ylim = c(minimo, 
                    massimo), xlab = "", ylab = "", col = "forestgreen", 
                    lwd = 4)
                if (cycle[indice] == "affine") 
                  points(sim, type = "b", ylim = c(minimo, massimo), 
                    xlab = "", ylab = "", col = "orange", lwd = 4)
            }
        }
        etichette <- rep(0, length(n.clust))
        for (i in 1:length(n.clust)) {
            etichette[i] <- paste("k =", n.clust[i])
        }
        box()
        boh <- 11
        axis(2)
        axis(1, at = 1:length(n.clust), labels = etichette, las = 0)
        ba <- c("Mean similarity indexes VS Num.clusters", paste("(center.method: ", 
            center.method, ")"))
        title(ba)
        leggenda <- cycle
        if (cycle[1] == "NOalignment") 
            leggenda[1] <- "without alignment"
        if (similarity.method.compare == "d0.L2" || similarity.method.compare == 
            "d1.L2" || similarity.method.compare == "d0.L2.centered" || 
            similarity.method.compare == "d1.L2.centered") {
            legend("topright", inset = 0.01, legend = leggenda, 
                col = colors, lwd = c(3, 3), cex = 0.6)
        }
        else {
            legend("bottomright", inset = 0.01, legend = leggenda, 
                col = colors, lwd = c(3, 3), cex = 0.6)
        }
    }
    if (plot.graph == 1 && length(n.clust) > 1) {
        toppo <- max(sim_NOalignment[1], sim_ONLYSHI[1], sim_ONLYDIL[1], 
            sim[1], na.rm = TRUE)
        uppo <- min(sim_NOalignment[1], sim_ONLYSHI[1], sim_ONLYDIL[1], 
            sim[1], na.rm = TRUE)
        lines(c(1, 1), c(uppo, toppo), col = "darkgrey", type = "l", 
            lty = "dashed", cex = 3, lwd = 3)
    }
    return(list(Result.NOalignment = Result.NOalignment, Result.shift = Result.shift, 
        Result.dilation = Result.dilation, Result.affine = Result.affine, 
        n.clust = n.clust, mean.similarity.NOalignment = sim_NOalignment, 
        mean.similarity.shift = sim_ONLYSHI, mean.similarity.dilation = sim_ONLYDIL, 
        mean.similarity.affine = sim))
}
