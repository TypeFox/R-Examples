# Estimate the active and non-active voxels based on the highest posterior density (HPD)
# of the coefficients simulated by the multilevel method.
# Plot the histogram of the posterior distribution of regression coefficient "vreg"}

regpostsim <-
function (pmeans, vreg, plot = T) 
{
    # HPD as in boa package
    boa.hpd <- function(x, alpha) {
        n <- length(x)
        m <- max(1, ceiling(alpha * n))
        y <- sort(x)
        a <- y[1:m]
        b <- y[(n - m + 1):n]
        i <- order(b - a)[1]
        structure(c(a[i], b[i]), names = c("Lower Bound", "Upper Bound"))
    }
    # Using Regression for Variable vreg
    # plotting examples;  posterior mean of var vreg
    pm1 <- pmeans[, 1]
    pm2 <- pmeans[, vreg]
    # cat("pm2 range\t", range(pm2), "\n")
    pm2.dens <- density(pm2)
    #---------------------------
    # HPD  based on posterior prob.
    hpd <- boa.hpd(pm2, 0.05)
    spmn <- which(pm2 < hpd[2])
    spma <- which(pm2 >= hpd[2])
    # cat("\nactive , non-active", length(spma), length(spmn), "\n")
    pm2.act <- pm2[spma]
    pm2.nact <- pm2[spmn]
    h1 <- 0.5 * max(pm2.dens$y)
    if (plot) {
        # x11(width = 7, height = 7)
        hist(pm2, prob = T, breaks = 15, col = "lightgray", xlab = paste("beta-", 
            vreg, sep = ""), main="Histogram of posterior mean values")
        lines(pm2.dens)
        print(pm2.dens)
        print("active range:")
        if (length(spma)) {
            print(range(pm2.act))
        }
        else {
            print("0")
        }
        if (length(spmn)) {
            print("non-active range:")
            print(range(pm2.nact))
        }
        else {
            print("0")
        }
        lines(c(hpd[1], hpd[1]), c(0, h1), lty = "dotted")
        lines(c(hpd[2], hpd[2]), c(0, h1), lty = "dotted")
        legend("topright", c("hpd interval (5%,95%)"), lty = c("dotted"))
        cat("hpd (95%)=      \t", hpd, "\n")
    }
    invisible(list(spma = spma, spmn = spmn))
}
