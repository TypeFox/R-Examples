show.estimate.hz <-
function (coordinates, path.mcmc.adm, burnin, angle = 0, plot.distruct = TRUE, 
    plot.mcmc = TRUE) 
{
    nindiv <- nrow(coordinates)
    param.adm <- read.table(file = paste(path.mcmc.adm, "parameters.hz.txt", 
        sep = ""))
    nit <- as.numeric(param.adm[param.adm[, 1] == "nit", 3])
    thinning <- as.numeric(param.adm[param.adm[, 1] == "thinning", 
        3])
    npop <- as.numeric(param.adm[param.adm[, 1] == "npop", 3])
    q.tmp <- as.matrix(read.table(paste(path.mcmc.adm, "q.txt", 
        sep = "")))
    q.mcmc <- array(dim = c(nindiv, npop, nit/thinning))
    iline <- 1
    for (iit in 1:(nit/thinning)) {
        for (iindiv in 1:nindiv) {
            q.mcmc[iindiv, , iit] <- q.tmp[iline, ]
            iline <- iline + 1
        }
    }
    q.est <- matrix(nrow = nindiv, ncol = npop)
    for (iindiv in 1:nindiv) {
        q.est[iindiv, ] <- apply(q.mcmc[iindiv, , -(1:burnin)], 
            1, mean)
    }
    if (plot.distruct) {
        nindiv <- nrow(coordinates)
        coord.proj <- numeric(nindiv)
        u <- matrix(c(cos(angle), sin(angle)), ncol = 1)
        P <- u %*% t(u)
        A <- matrix(c(cos(-angle), -sin(-angle), sin(-angle), 
            cos(-angle)), byrow = TRUE, ncol = 2)
        coord.proj <- t(A %*% (P %*% t(coordinates)))[, 1]
        q.sort <- q.est[order(coord.proj), ]
        palette(terrain.colors(npop))
        plot((1:nindiv)/nindiv, (1:nindiv)/nindiv, type = "n", 
            xlab = "Spatial coordinate along a one-dimensional axis", 
            ylab = "Admixture proportions", ylim = c(0, 1), lab = c(1, 
                10, 7))
        for (iindiv in 1:(nindiv - 1)) {
            for (ipop in 1:(npop - 1)) {
                polygon(c(iindiv/nindiv, (iindiv + 1)/nindiv, 
                  (iindiv + 1)/nindiv, iindiv/nindiv), c(sum(q.sort[iindiv, 
                  1:ipop]), sum(q.sort[iindiv, 1:ipop]), sum(q.sort[iindiv, 
                  1:(ipop + 1)]), sum(q.sort[iindiv, 1:(ipop + 
                  1)])), col = ipop + 1, border = NA)
            }
            polygon(c(iindiv/nindiv, (iindiv + 1)/nindiv, (iindiv + 
                1)/nindiv, iindiv/nindiv), c(0, 0, q.sort[iindiv, 
                1], q.sort[iindiv, 1]), col = 1, border = NA)
        }
        palette("default")
    }
    if (plot.mcmc) {
        a.mcmc <- as.matrix(read.table(paste(path.mcmc.adm, "a.txt", 
            sep = "")))
        b.mcmc <- as.matrix(read.table(paste(path.mcmc.adm, "b.txt", 
            sep = "")))
        c.mcmc <- as.matrix(read.table(paste(path.mcmc.adm, "c.txt", 
            sep = "")))
        a.est <- apply(a.mcmc[-(1:burnin), ], 2, mean)
        b.est <- apply(b.mcmc[-(1:burnin), ], 2, mean)
        c.est <- apply(c.mcmc[-(1:burnin), ], 2, mean)
        dev.new()
        par(mfrow = c(1, 3))
        plot(a.mcmc[, 1], type = "n", xlab = "Index of saved iteration", 
            ylab = "a")
        for (ipop in 1:npop) {
            lines(a.mcmc[, ipop], lty = ipop, col = ipop)
            abline(h = a.est[ipop], lty = ipop, col = ipop)
        }
        plot(b.mcmc[, 1], type = "n", xlab = "Index of saved iteration", 
            ylab = "b")
        for (ipop in 1:npop) {
            lines(b.mcmc[, ipop], lty = ipop, col = ipop)
            abline(h = b.est[ipop], lty = ipop, col = ipop)
        }
        plot(c.mcmc[, 1], type = "n")
        for (ipop in 1:npop) {
            lines(c.mcmc[, ipop], lty = ipop, col = ipop)
            abline(h = c.est[ipop], lty = ipop, col = ipop)
        }
    }
}
