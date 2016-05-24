loglik.norm.plot <- function (X, parameter = c("mu", "sigma.sq"), poss = NULL, plot.likfunc = TRUE, 
    plot.density = TRUE, plot.calc = FALSE, xlab = NULL, ylab = NULL, 
    conv = 0.01, anim = TRUE, est.col = 2, density.leg = TRUE, cex.leg = 0.9, interval = 0.01, ...) 
{
    possibilities <- poss
    Var.MLE <- function(e.i) sum((e.i - mean(e.i))^2)/length(e.i)
    if (is.null(poss)) {
        if (parameter == "mu") 
            possibilities <- seq(mean(X) - 2 * sd(X), mean(X) + 
                2 * sd(X), conv)
        if (parameter == "sigma.sq") 
            possibilities <- seq(0.5 * Var.MLE(X), 2 * Var.MLE(X), 
                conv)
        logl <- seq(min(possibilities), max(possibilities), conv)
    }
    if (!is.null(poss)) 
        logl <- poss
    for (i in 1:length(possibilities)) {
        if (parameter == "mu") {
            sdM <- sqrt(Var.MLE(X))
            logl[i] <- log(prod(dnorm(X, mean = possibilities[i], 
                sd = sdM)))
        }
        if (parameter == "sigma.sq") {
            mu <- mean(X)
            logl[i] <- log(prod(dnorm(X, mean = mu, sd = sqrt(possibilities[i]))))
        }
    }
    max.p <- ifelse(is.null(poss), ifelse(parameter == "mu", 
        mean(X), Var.MLE(X)), possibilities[logl == max(logl)][1])
    MLE <- ifelse(parameter == "mu", mean(X), Var.MLE(X))
    Max.lik <- log(prod(dnorm(X, mean = mean(X), sd = sqrt(Var.MLE(X)))))
    if (anim == FALSE) {
        if (plot.likfunc == TRUE & plot.density == TRUE) {
            dev.new(height = 4, width = 8)
            par(mfrow = c(1, 2), mar = c(4.4, 4.5, 1, 0.5), cex = 0.9)
            layout(matrix(c(1, 2), 1, 2, byrow = TRUE))
        }
        if (plot.likfunc == TRUE) {
            plot(possibilities, logl, type = "l", ylab = ifelse(is.null(ylab), 
                "Normal log-likelihood function", ylab), xlab = ifelse(is.null(xlab), 
                ifelse(parameter == "mu", expression(paste("Estimates for ", 
                  mu)), expression(paste("Estimates for ", sigma^2))), 
                xlab), ...)
            legend("topright", legend = bquote(paste("ML est. = ", 
                .(MLE))), cex = cex.leg, bty = "n")
            abline(v = MLE, lty = 2, col = est.col)
        }
        if (plot.density == TRUE) {
            x <- NULL
            rm(x)
            curve(dnorm(x, mean = mean(X), sd = sqrt(Var.MLE(X))), 
                from = min(X), to = max(X), xlab = expression(italic(x)), 
                ylab = expression(paste(italic(f), "(", italic(x), 
                  ")", sep = "")))
              legend("topright",legend=c(paste("X~N(",bquote(.(round(mean(X),1))),", ",bquote(.(round(Var.MLE(X),1))),")",sep=""),
              paste("Loglik = ",bquote(.(round(Max.lik,2))))),bty="n",cex=cex.leg)
                
             
            if(density.leg == TRUE){legend("topleft", col = 2, lty = 2, legend = "Obs. density", 
                bty = "n", cex = cex.leg)}
            segments(x0 = X, y0 = rep(0, length(X)), x1 = X, 
                y1 = dnorm(X, mean = mean(X), sd = sqrt(Var.MLE(X))), 
                col = est.col, lty = 2)
        }
    }
    if (anim == TRUE) {
        nm <- which(logl == max(logl))[1]
        if (plot.likfunc == TRUE & plot.density == TRUE) 
            dev.new(height = 4, width = 8)
        par(mfrow = c(1, 2), mar = c(4.4, 4.5, 1, 0.5), cex = 0.9)
        for (i in 1:(nm - 1)) {
            dev.hold()
            if (plot.likfunc == TRUE) {
                plot(possibilities, logl, type = "n", ylab = ifelse(is.null(ylab), 
                  "Normal log-likelihood function", ylab), xlab = ifelse(is.null(xlab), 
                  ifelse(parameter == "mu", expression(paste("Estimates for ", 
                    mu)), expression(paste("Estimates for ", 
                    sigma^2))), xlab), ...)
                arrows(possibilities[i], logl[i], possibilities[i + 
                  1], logl[i + 1], col = est.col, length = 0.15, lwd = 1)
                points(possibilities[1:i], logl[1:i], lty = 2, 
                  col = est.col, lwd = 1, type = "l")
                if (i == (nm - 1)) {
                  points(possibilities, logl, type = "l")
                  abline(v = MLE, lty = 2, col = est.col)
                  legend("topright", legend = bquote(paste("ML est. = ", 
                    .(MLE))), cex = cex.leg, bty = "n")
                }
            }
            if (plot.density == TRUE) {
                x <- NULL
                rm(x)
                curve(dnorm(x, mean = ifelse(parameter == "mu", 
                  possibilities[i], mu), sd = ifelse(parameter == 
                  "mu", sdM, possibilities[i])), from = min(X), 
                  to = max(X), xlab = expression(italic(x)), 
                  ylab = expression(paste(italic(f), "(", italic(x), 
                    ")", sep = "")))
                segments(x0 = X, y0 = rep(0, length(X)), x1 = X, 
                  y1 = dnorm(X, mean = ifelse(parameter == "mu", 
                    possibilities[i], mu), sd = ifelse(parameter == 
                    "mu", sdM, possibilities[i])), col = est.col, lty = 2)
                if (i != (nm - 1)) {
                  legend("topright", legend = c(ifelse(parameter == 
                    "mu", as.expression(bquote(paste(mu, " = ", 
                    .(round(possibilities[i], 2))))), as.expression(bquote(paste(sigma^2, 
                    " = ", .(round(possibilities[i], 2)))))), 
                    as.expression(bquote(paste("Loglik = ", .(round(logl[i], 
                      2)))))), bty = "n", cex = cex.leg)
                }
                if (i == (nm - 1)) {
                  legend("topright", legend = c(ifelse(parameter == 
                    "mu", as.expression(bquote(paste(mu, " = ", 
                    .(MLE)))), as.expression(bquote(paste(sigma^2, 
                    " = ", .(MLE))))), as.expression(bquote(paste("Loglik = ", 
                    .(round(Max.lik, 2)))))), bty = "n", cex = cex.leg)
                }
            }
            dev.flush()
            Sys.sleep(interval)
        }
        if (plot.calc == TRUE) {
            dev.new(width = 4, height = 4, ypos = 225)
            par(mar = c(4.4, 4.5, 1, 0.5), cex = 0.9)
            X <- sort(X)
            c <- ifelse(length(X) < 200, round(200/length(X), 
                0), 1)
            x0 = rep(X, each = c)
            y0 = rep(rep(0, length(X)), each = c)
            y1 = rep(dnorm(X, mean = ifelse(parameter == "mu", 
                max.p, mu), sd = ifelse(parameter == "mu", sdM, 
                max.p)), each = c)
            loglik <- log(dnorm(X, mean = ifelse(parameter == 
                "mu", max.p, mu), sd = ifelse(parameter == "mu", 
                sdM, max.p)))
            cumlik <- rep(cumsum(loglik), each = c)
            for (i in 1:(length(X) * c)) {
                dev.hold()
                x <- NULL
                rm(x)
                curve(dnorm(x, mean = ifelse(parameter == "mu", 
                  max.p, mu), sd = ifelse(parameter == "mu", 
                  sdM, max.p)), from = min(X), to = max(X), xlab = expression(italic(x)), 
                  ylab = expression(paste(italic(f), "(", italic(x), 
                    ")", sep = "")))
                if(density.leg == TRUE){legend("topleft", col = 2, lty = 2, legend = "Obs. density", 
                  bty = "n", cex = 0.9)}
                legend("topright", legend = c(paste("X~N(", bquote(.(round(mean(X), 
                  1))), ", ", bquote(.(round(Var.MLE(X), 1))), 
                  ")", sep = ""), paste("Loglik = ", bquote(.(round(cumlik[i], 
                  2))))), bty = "n", cex = cex.leg)
                segments(x0[1:i], y0[1:i], x0[1:i], y1[1:i], 
                  col = 2, lty = 2)
                dev.flush()
                Sys.sleep(interval)
            }
        }
    }
}
