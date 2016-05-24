fittedPlot <- function (object, ..., x = NULL, color = TRUE, line.type = FALSE, xlab=NULL) 
{
    if (length(list(...))) {
        object <- list(object, ...)
        nobj <- length(object)
        isgamlss <- unlist(lapply(object, is.gamlss))
        if (!any(isgamlss)) 
            stop("some of the objects are not gamlss")
        nopar <- as.numeric(lapply(object, function(x) length(x$parameters)))
        param <- lapply(object, function(x) x$parameters)
        xvar <- if (is.null(xlab)) deparse(substitute(x)) else xlab
        x.o <- x[order(x)]
        let <- c("(a)", "(b)", "(c)", "(d)")
        index.mu <- index.sigma <- index.nu <- index.tau <- 1
        for (ii in 1:nobj) {
            if ("mu" %in% param[[ii]]) {
                mu.o <- fitted(object[[ii]], "mu")[order(x)]
                if (all(abs(mu.o - mu.o[1]) < 1e-04)) 
                  mu.o <- rep(mu.o[1], length(mu.o))
                mu.mat <- if (index.mu == 1) 
                  mu.o
                else cbind(mu.mat, mu.o)
                index.mu <- 2
            }
            if ("sigma" %in% param[[ii]]) {
                sigma.o <- fitted(object[[ii]], "sigma")[order(x)]
                if (all(abs(sigma.o - sigma.o[1]) < 1e-04)) 
                  sigma.o <- rep(sigma.o[1], length(sigma.o))
                sigma.mat <- if (index.sigma == 1) 
                  sigma.o
                else cbind(sigma.mat, sigma.o)
                index.sigma <- 2
            }
            if ("nu" %in% param[[ii]]) {
                nu.o <- fitted(object[[ii]], "nu")[order(x)]
                if (all(abs(nu.o - nu.o[1]) < 1e-04)) 
                  nu.o <- rep(nu.o[1], length(nu.o))
                nu.mat <- if (index.nu == 1) 
                  nu.o
                else cbind(nu.mat, nu.o)
                index.nu <- 2
            }
            if ("tau" %in% param[[ii]]) {
                tau.o <- fitted(object[[ii]], "tau")[order(x)]
                if (all(abs(tau.o - tau.o[1]) < 1e-04)) 
                  tau.o <- rep(tau.o[1], length(tau.o))
                tau.mat <- if (index.tau == 1) 
                  tau.o
                else cbind(tau.mat, tau.o)
                index.tau <- 2
            }
        }
        col <- 1
        vvv <- switch(max(nopar), c(1, 1), c(2, 1), c(3, 1), 
            c(2, 2))
        op <- par(mfrow = vvv, mar = par("mar") + c(0, 1, 0, 
            0), col.axis = "blue4", col.main = "blue4", col.lab = "blue4", 
            cex = 0.5, cex.lab = 1.3, cex.axis = 1, cex.main = 1.3)
        if ("mu" %in% unlist(param)) {
            ncolmu <- if (is.null(dim(mu.mat))) 
                1
            else dim(mu.mat)[2]
            if (ncolmu == 1) 
                mu.mat <- matrix(mu.mat, ncol = 1)
            if (color == TRUE) 
                col <- 3
            ltype <- 1
            plot(x.o, mu.mat[, 1], xlab = xvar, ylab = "mu", 
                main = let[1], col = "darkgreen", col.axis = "mediumblue", 
                type = "n", ylim = range(mu.mat), frame.plot = TRUE)
            for (ii in 1:ncolmu) {
                lines(x.o, mu.mat[, ii], col = col, lty = ltype)
                if (color == TRUE) 
                  col <- col + 1
                if (line.type == TRUE) 
                  ltype <- ltype + 1
            }
        }
        if ("sigma" %in% unlist(param)) {
            ncolsigma <- if (is.null(dim(sigma.mat))) 
                1
            else dim(sigma.mat)[2]
            if (ncolsigma == 1) 
                sigma.mat <- matrix(sigma.mat, ncol = 1)
            if (color == TRUE) 
                col <- 3
            ltype <- 1
            plot(x.o, sigma.mat[, 1], xlab = xvar, ylab = "sigma", 
                main = let[2], col = "darkgreen", col.axis = "mediumblue", 
                type = "n", ylim = range(sigma.mat), frame.plot = TRUE)
            for (ii in 1:ncolsigma) {
                lines(x.o, sigma.mat[, ii], col = col, lty = ltype)
                if (color == TRUE) 
                  col <- col + 1
                if (line.type == TRUE) 
                  ltype <- ltype + 1
            }
        }
        if ("nu" %in% unlist(param)) {
            ncolnu <- if (is.null(dim(nu.mat))) 
                1
            else dim(nu.mat)[2]
            if (ncolnu == 1) 
                nu.mat <- matrix(nu.mat, ncol = 1)
            if (color == TRUE) 
                col <- 3
            ltype <- 1
            plot(x.o, nu.mat[, 1], xlab = xvar, ylab = "nu", 
                main = let[3], col = "darkgreen", col.axis = "mediumblue", 
                type = "n", ylim = range(nu.mat), frame.plot = TRUE)
            for (ii in 1:ncolnu) {
                lines(x.o, nu.mat[, ii], col = col, lty = ltype)
                if (color == TRUE) 
                  col <- col + 1
                if (line.type == TRUE) 
                  ltype <- ltype + 1
            }
        }
        if ("tau" %in% unlist(param)) {
            ncoltau <- if (is.null(dim(tau.mat))) 
                1
            else dim(tau.mat)[2]
            if (ncoltau == 1) 
                tau.mat <- matrix(tau.mat, ncol = 1)
            if (color == TRUE) 
                col <- 3
            ltype <- 1
            plot(x.o, tau.mat[, 1], xlab = xvar, ylab = "tau", 
                main = let[4], col = "darkgreen", col.axis = "mediumblue", 
                type = "n", ylim = range(tau.mat), frame.plot = TRUE)
            for (ii in 1:ncoltau) {
                lines(x.o, tau.mat[, ii], col = col, lty = ltype)
                if (color == TRUE) 
                  col <- col + 1
                if (line.type == TRUE) 
                  ltype <- ltype + 1
            }
        }
        par(op)
    }
    else {
        if (!is.gamlss(object)) 
            stop(paste("This is not an gamlss object", "\n", 
                ""))
        if (is.null(x)) 
            stop(paste("The x-variable argument is not specified", 
                "\n", ""))
        nopar <- length(object$parameters)
        param <- object$parameters
        xvar <- if (is.null(xlab)) deparse(substitute(x)) else xlab
        x.o <- x[order(x)]
        let <- c("(a)", "(b)", "(c)", "(d)")
        if ("mu" %in% object$parameters) {
            vvv <- c(1, 1)
            mu.o <- fitted(object, "mu")[order(x)]
            if (all(abs(mu.o - mu.o[1]) < 1e-04)) 
                mu.o <- rep(mu.o[1], length(mu.o))
            mat <- cbind(mu.o)
        }
        if ("sigma" %in% object$parameters) {
            vvv <- c(2, 1)
            sigma.o <- fitted(object, "sigma")[order(x)]
            if (all(abs(sigma.o - sigma.o[1]) < 1e-04)) 
                sigma.o <- rep(sigma.o[1], length(sigma.o))
            mat <- cbind(mat, sigma.o)
        }
        if ("nu" %in% object$parameters) {
            vvv <- c(3, 1)
            nu.o <- fitted(object, "nu")[order(x)]
            if (all(abs(nu.o - nu.o[1]) < 1e-04)) 
                nu.o <- rep(nu.o[1], length(nu.o))
            mat <- cbind(mat, nu.o)
        }
        if ("tau" %in% object$parameters) {
            vvv <- c(2, 2)
            tau.o <- fitted(object, "tau")[order(x)]
            if (all(abs(tau.o - tau.o[1]) < 1e-04)) 
                tau.o <- rep(tau.o[1], length(tau.o))
            mat <- cbind(mat, tau.o)
        }
        op <- par(mfrow = vvv, mar = par("mar") + c(0, 1, 0, 
            0), col.axis = "blue4", col.main = "blue4", col.lab = "blue4", 
            cex = 0.5, cex.lab = 1.3, cex.axis = 1, cex.main = 1.3)
        for (ii in 1:nopar) {
            pp <- as.character(param[ii])
            if (color == TRUE) 
                col <- "darkgreen"
            else col <- "black"
            plot(x.o, mat[, ii], xlab = xvar, ylab = pp, main = let[ii], 
                col = col, col.axis = "mediumblue", type = "l", 
                frame.plot = TRUE)
        }
        par(op)
    }
}
