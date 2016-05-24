maxDiffCDF <- function(res1, res2, which = c("MLE", "smooth"), n.grid = 500){

# calculate maximal difference between
# F_x and F_y via zeros of phi_x - phi_y
x <- res1$knots
y <- res2$knots

z <- c(x, y)
d <- diff(range(z))

test.stat <- rep(NA, 2)
names(test.stat) <- c("log-concave", "smooth log-concave")
loc <- rep(NA, 2)
names(loc) <- names(test.stat)


if ("MLE" %in% which){
    ## for log-concave density CDF
    l <- max(min(x), min(y))
    r <- min(max(x), max(y))
    z <- sort(unique(c(z[(z >= l) & (z <= r)])))
    N <- length(z)
    if (r <= l){
        test.stat <- 1
        loc <- r
        } else {

    mat1 <- matrix(NA, nrow = N, ncol = 3)
    mat2 <- mat1

    # compute values of functions at z
    for (i in 1:N){
        mat1[i, ] <- evaluateLogConDens(z[i], res1, which = 1)[, "log-density"]
        mat2[i, ] <- evaluateLogConDens(z[i], res2, which = 1)[, "log-density"]
    }

    # find places of crossings
    x0 <- matrix(c(z[1], abs(evaluateLogConDens(z[1], res1, which = 3)[, "CDF"] - evaluateLogConDens(z[1], res2, which = 3)[, "CDF"] )), ncol = 2)   # difference at z[1]
    i <- 1:(length(z) - 1)
    i <- i[sign(mat1[i, 1] - mat2[i, 1]) != sign(mat1[i + 1, 1] - mat2[i + 1, 1])]

    if (length(i) > 0){

        # calculate exact places of crossings by setting
        # phi_x(x) = phi_y(x) and solving for x
        zk <- z[i]
        zk1 <- z[i + 1]
        phi1k <- mat1[i, 1]
        phi1k1 <- mat1[i + 1, 1]
        phi2k <- mat2[i, 1]
        phi2k1 <- mat2[i + 1, 1]
        s1 <- (phi1k1 - phi1k) / (zk1 - zk)
        s2 <- (phi2k1 - phi2k) / (zk1 - zk)

        tmp <- sort(unique(c(zk[s1 == s2], (zk + (phi2k - phi1k) / (s1 - s2))[s1 != s2])))
        for (j in 1:length(tmp)){x0 <- rbind(x0, c(tmp[j], abs(evaluateLogConDens(tmp[j], res1, which = 3)[, "CDF"] - evaluateLogConDens(tmp[j], res2, which = 3)[, "CDF"])))}
    }

    # add difference of distribution functions at z[N]
    x0 <- rbind(x0, c(z[N], abs(evaluateLogConDens(z[N], res1, which = 3)[, "CDF"] - evaluateLogConDens(z[N], res2, which = 3)[, "CDF"])))

    test.stat[1] <- max(x0[, 2])
    loc[1] <- x0[x0[, 2] == test.stat[1], 1]
    }
}


if ("smooth" %in% which){
    ## for smooth log-concave density CDF
    xs <- seq(min(x, y) - 0.1 * d, max(x, y) + 0.1 * d, length.out = n.grid)
    F1smooth <- rep(NA, length(xs)); F2smooth <- F1smooth
    for (i in 1:length(xs)){
        F1smooth[i] <- evaluateLogConDens(xs[i], res1, which = 4)[, "smooth.density"] 
        F2smooth[i] <- evaluateLogConDens(xs[i], res2, which = 4)[, "smooth.density"] 
    }
    diff <- F1smooth - F2smooth
    ind1 <- abs(diff(sign(diff)))

    left <- xs[c(ind1, 0) == 2]
    right <- xs[c(0, ind1) == 2]

    locs <- rep(NA, length(left))
    test.stats <- locs
    zero.fun <- function(q, res1, res2){return(evaluateLogConDens(q, res1, which = 4)[, "smooth.density"]  - evaluateLogConDens(q, res2, which = 4)[, "smooth.density"] )}

    for (i in 1:length(locs)){
        locs[i] <- uniroot(f = zero.fun, interval = c(left[i], right[i]), res1, res2)$root
        test.stats[i] <- abs(evaluateLogConDens(locs[i], res1, which = 5)[, "smooth.CDF"] - evaluateLogConDens(locs[i], res2, which = 5)[, "smooth.CDF"])
    }

    loc[2] <- locs[test.stats == max(test.stats)]
    test.stat[2] <- test.stats[test.stats == max(test.stats)]
}

return(list("test.stat" = as.numeric(test.stat), "location" = as.numeric(loc)))
}







