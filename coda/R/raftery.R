"raftery.diag" <-
function (data, q = 0.025, r = 0.005, s = 0.95, converge.eps = 0.001) 
{
  if (is.mcmc.list(data)) 
        return(lapply(data, raftery.diag, q, r, s, converge.eps))
    data <- as.mcmc(data)
    resmatrix <- matrix(nrow = nvar(data), ncol = 4, dimnames = list(varnames(data, 
        allow.null = TRUE), c("M", "N", "Nmin", "I")))
    phi <- qnorm(0.5 * (1 + s))
    nmin <- as.integer(ceiling((q * (1 - q) * phi^2)/r^2))
    if (nmin > niter(data)) 
        resmatrix <- c("Error", nmin)
    else for (i in 1:nvar(data)) {
        #          First need to find the thinning parameter kthin 
        # 
        if (is.matrix(data)) {
            quant <- quantile(data[, i, drop = TRUE], probs = q)
            dichot <- mcmc(data[, i, drop = TRUE] <= quant, start = start(data), 
                end = end(data), thin = thin(data))
        }
        else {
            quant <- quantile(data, probs = q)
            dichot <- mcmc(data <= quant, start = start(data), 
                end = end(data), thin = thin(data))
        }
        kthin <- 0
        bic <- 1
        while (bic >= 0) {
            kthin <- kthin + thin(data)
            testres <- as.vector(window.mcmc(dichot, thin = kthin))
            testres <- factor(testres, levels=c(FALSE,TRUE))
            newdim <- length(testres)
            testtran <- table(testres[1:(newdim - 2)], testres[2:(newdim - 
                1)], testres[3:newdim])
            testtran <- array(as.double(testtran), dim = dim(testtran))
            g2 <- 0
            for (i1 in 1:2) {
                for (i2 in 1:2) {
                  for (i3 in 1:2) {
                    if (testtran[i1, i2, i3] != 0) {
                      fitted <- (sum(testtran[i1, i2, 1:2]) * 
                        sum(testtran[1:2, i2, i3]))/(sum(testtran[1:2, 
                        i2, 1:2]))
                      g2 <- g2 + testtran[i1, i2, i3] * log(testtran[i1, 
                        i2, i3]/fitted) * 2
                    }
                  }
                }
            }
            bic <- g2 - log(newdim - 2) * 2
        }
        #
        # then need to find length of burn-in and No of iterations for required precision 
        # 
        finaltran <- table(testres[1:(newdim - 1)], testres[2:newdim])
        alpha <- finaltran[1, 2]/(finaltran[1, 1] + finaltran[1, 
            2])
        beta <- finaltran[2, 1]/(finaltran[2, 1] + finaltran[2, 
            2])
        tempburn <- log((converge.eps * (alpha + beta))/max(alpha, 
            beta))/(log(abs(1 - alpha - beta)))
        nburn <- as.integer(ceiling(tempburn) * kthin)
        tempprec <- ((2 - alpha - beta) * alpha * beta * phi^2)/(((alpha + 
            beta)^3) * r^2)
        nkeep <- as.integer(ceiling(tempprec) * kthin)
        iratio <- (nburn + nkeep)/nmin
        resmatrix[i, 1] <- nburn
        resmatrix[i, 2] <- nkeep + nburn
        resmatrix[i, 3] <- nmin
        resmatrix[i, 4] <- signif(iratio, digits = 3)
    }
    y <- list(params = c(r = r, s = s, q = q), resmatrix = resmatrix)
    class(y) <- "raftery.diag"
    return(y)
}
"print.raftery.diag" <-
function (x, digits = 3, ...) 
{
    cat("\nQuantile (q) =", x$params["q"])
    cat("\nAccuracy (r) = +/-", x$params["r"])
    cat("\nProbability (s) =", x$params["s"], "\n")
    if (x$resmatrix[1] == "Error") 
        cat("\nYou need a sample size of at least", x$resmatrix[2], 
            "with these values of q, r and s\n")
    else {
        out <- x$resmatrix
        for (i in ncol(out)) out[, i] <- format(out[, i], digits = digits)
        out <- rbind(matrix(c("Burn-in ", "Total", "Lower bound ", 
            "Dependence", "(M)", "(N)", "(Nmin)", "factor (I)"), 
            byrow = TRUE, nrow = 2), out)
        if (!is.null(rownames(x$resmatrix))) 
            out <- cbind(c("", "", rownames(x$resmatrix)), out)
        dimnames(out) <- list(rep("", nrow(out)), rep("", ncol(out)))
        print.default(out, quote = FALSE, ...)
        cat("\n")
    }
    invisible(x)
}
