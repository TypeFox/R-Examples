library(robustbase)

n <- 1:50
(qnn <- sapply(n, function(n)Qn(1:n, const=1)))
plot(n, qnn, type = 'b', col = 2,
     ylab = "Qn", main = "Qn(1:n) [unscaled]")

(snn <- sapply(n, function(n)Sn(1:n, const=1)))
plot(n, snn, type = 'b', col = 2,
     ylab = "Sn", main = "Sn(1:n) [unscaled]")

matplot(n, cbind(qnn, snn),type = 'b',
        ylab = "Qn & Sn", main = "Qn(1:n) & Sn(1:n) [unscaled]")

(sdn <- c(1, sapply(n[-1], function(n)sd(1:n)/n)))
## sd(1) => NA

for(Sample in c(function(n) ppoints(n),
                function(n) qnorm(ppoints(n))))
{
    ##mult.fig(2) :
    op <- par(mfrow=c(2,1), mgp = c(1.5,.6,0), mar = .1 + c(4,4,2,1))
    for(N in c(50, 200)) {
        n <- 1:N
        sdn <- c(1, sapply(n[-1], function(m)sd(Sample(m))))
        r <- cbind(Qn = sapply(n, function(m)Qn(Sample(m))),
                   Sn = sapply(n, function(m)Sn(Sample(m)))) / sdn
        matplot(n, r, type = 'b', col = 2:3, lty = 1, ylab = "Qn & Sn",
                main = "Qn(Sample(n)) & Sn(..) [consistently scaled]")
        legend(.85*N, 0.4, c("Qn()", "Sn()"), col = 2:3, lty = 1,
               pch = paste(1:2))
        abline(h=1, col = "gray", lty = 2)
    }
    par(op)

    ## Hmm,  the above does not look 100% consistent to *my* eyes...
    ## Investigate:

    matplot(n, r, ylim = c(0.9, 1.1), type = 'b', col = 2:3, lty = 1)
    abline(h=1, col = "gray", lty = 2)

    matplot(n, r^2, ylim = c(0.7, 1.3), type = 'b', col = 2:3, lty = 1)
    abline(h=1, col = "gray", lty = 2)
}
rownames(r) <- paste(n)
r
