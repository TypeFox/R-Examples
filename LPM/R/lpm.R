lpm <-
function (x, p, q, n, smean, svar, outer = 0, prob = 0.95, fre = 365, 
    fractional = F, Plag = 20, lsign = 0.05, n1 = 399, trasfo = F, 
    des = T, rain = F, graph = F) 
{if ((n == 0)) rain = FALSE
    lsign = 1 - lsign
    if (fractional) 
        qint = qt((1 - (1 - prob)/2), (length(x) - (p + q + 1) - 
            1))
    else qint = qt((1 - (1 - prob)/2), (length(x) - (p + q) - 
        1))
    if (n > 0) 
        simul = T
    else simul = F
    options(warn = -1)
    if (trasfo) {
        x <- log(x)
    }
    if (des) {
        x1 <- Stagiona(x, x, fre, outer, smean, svar, 0, par)
        if (graph) {
            dev.new()
            plot(x1$M, type = "l", main = "Mean", xlab = "", 
                ylab = "")
            dev.new()
            plot(x1$V, type = "l", main = "St. Deviation", xlab = "", 
                ylab = "")
        }
        x1 = x1$des
    }
    else x1 = x
    mx = mean(x1)
    x1 = x1 - mx
    if (fractional == F) {
        x2 <- modarma(x1, p, q, Plag, n1, n, simul, graph = graph, 
            lsign = lsign)
        cat("Portmonteau Test :", x2$BIC$Q, x2$BIC$DistrChi, 
            "\n")
        cat("BIC:", x2$BIC$BIC, "\n")
        if (p > 0) 
            cat("Parameters: AR", x2$para$ar, "\n")
        cat("MA", x2$para$ma, "\n")
        cat(prob * 100, "% Confidence Bounds", "\n")
        if (x2$para$ar != 0) {
            for (b in 1:length(x2$para$ar)) cat("[", x2$para$ar[b] - 
                qint * x2$stdev[b], ";", x2$para$ar[b] + qint * 
                x2$stdev[b], "]", "\n")
            if (x2$para$ma != 0) 
                for (c in 1:length(x2$para$ma)) cat("[", x2$para$ma[c] - 
                  qint * x2$stdev[(length(x2$para$ar) + c)], 
                  ";", x2$para$ma[c] + qint * x2$stdev[(length(x2$para$ar) + 
                    c)], "]", "\n")
        }
        else for (c in 1:length(x2$para$ma)) cat("[", x2$para$ma[c] - 
            qint * x2$stdev[c], ";", x2$para$ma[c] + qint * x2$stdev[c], 
            "]", "\n")
        npar <- 0
    }
    if (fractional) {
        x2 <- modfarma(x1, p, q, Plag, n, n1, simul, graph = graph, 
            lsign = lsign)
        cat("Portmanteau Test :", x2$BIC$Q, x2$BIC$DistrChi, 
            "\n")
        cat("BIC:", x2$BIC$BIC, "\n")
        cat("Parameters: d", x2$parameter$d, "\n")
        cat("AR", x2$parameter$AR, "\n")
        cat("MA", x2$parameter$MA, "\n")
        cat(prob * 100, "% Confidence Bounds", "\n")
        cat("[", x2$parameter$d - qint * x2$stdev[1], ";", x2$parameter$d + 
            qint * x2$stdev[1], "]", "\n")
        if ((p + q) > 0) {
            if (x2$parameter$AR != 0) {
                for (b in 1:length(x2$parameter$AR)) cat("[", 
                  x2$parameter$AR[b] - qint * x2$stdev[(b + 1)], 
                  ";", x2$parameter$AR[b] + qint * x2$stdev[(b + 
                    1)], "]", "\n")
                if (x2$parameter$MA != 0) 
                  for (c in 1:length(x2$parameter$MA)) cat("[", 
                    x2$parameter$MA[c] - qint * x2$stdev[(length(x2$parameter$AR) + 
                      c + 1)], ";", x2$parameter$MA[c] + qint * 
                      x2$stdev[(length(x2$parameter$AR) + c + 
                        1)], "]", "\n")
            }
            else {
                for (c in 1:length(x2$parameter$MA)) cat("[", 
                  x2$parameter$MA[c] - qint * x2$stdev[(c + 1)], 
                  ";", x2$parameter$MA[c] + qint * x2$stdev[(c + 
                    1)], "]", "\n")
            }
        }
        npar <- 0
    }
    if (simul) {
        x2a <- x2$simulazioni
        xfin <- list()
        for (i in 1:n) {
            x2a[[i]] = x2a[[i]] + mx
            if (des) {
                x3 <- Stagiona(x, x2a[[i]], fre, outer, smean, 
                  svar, 1, npar)
                xfin[[i]] <- x3$stag
            }
            else xfin <- x2a
            if (trasfo) {
                xfin[[i]] <- exp(xfin[[i]])
            }
        }
    }
    if (rain) {
        for (j in 1:n) {
            cat("\nRain adaptor sim", j, "\n")
            xfin[[j]] = rain.adapt(as.matrix(x), xfin[[j]], 1)
        }
    }
    results <- list()
    results$res <- x2$residui
    results$para <- x2$para
    results$stdev <- x2$stdev
    results$Bic <- x2$BIC$BIC
    if (simul) {
        results$simdes = x2a
        results$sim <- xfin
    }
    return(results)
}
