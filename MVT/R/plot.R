envelope <- function(object, reps = 50, conf = 0.95, plot.it = TRUE)
{   # simulated envelope
    envel <- function(n, mean, Sigma, eta, reps, conf)
    {
        conf <- 1 - conf
        # initialize progress bar
        cat("  Progress:\n")
        pb <- txtProgressBar(min = 0, max = reps, style = 3)
        elims <- matrix(0, nrow = n, ncol = reps)
        for (i in 1:reps) {
            x <- rmt(n, mean = mean, Sigma = Sigma, eta = eta)
            fit <- studentFit(x, family = Student(eta = eta))
            z <- wilson.hilferty(fit, eta = fit$eta)
            elims[,i] <- sort(z)
            # update progress bar
            setTxtProgressBar(pb, i)
        }
        close(pb)
        band <- matrix(0, nrow = n, ncol = 2)
        for (i in 1:n)
            band[i,] <- quantile(elims[i,], probs = c(conf / 2, 1 - conf / 2))
        band
    }
    
    n <- object$dims[1]
    z <- wilson.hilferty(object, eta = object$eta)
    
    if (plot.it) {
        band  <- envel(n, object$center, object$Scatter, object$eta, reps, conf)
        ylim <- range(z, band)
        qqnorm(z, ylim = ylim, main = "Transformed distances Q-Q plot")
        par(new = TRUE)
        qqnorm(band[,1], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
        par(new = TRUE)
        qqnorm(band[,2], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
    }
    invisible(list(transformed = z, envelope = band))
}

wilson.hilferty <- function(x, center, cov, eta = 0)
{   # Wilson-Hilferty transformation
    if(missing(center) && missing(cov)) {
        if (is.null(x$distances) && is.null(x$dims))
            stop("x is not a valid object")
        distances <- x$distances
        p <- x$dims[2]
    }
    else {
        distances <- mahalanobis(x, center, cov)
        p <- nrow(cov)
    }
        
    y <- distances / (1 - 2 * eta)
    y <- y / p
    z <- ((1 - 2 * eta / 9) * y^(1/3) - (1 - 2 / (9 * p))) / sqrt(2 * eta * y^(2/3) / 9 + 2 / (9 * p))
    z
}
