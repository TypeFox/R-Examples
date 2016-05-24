## 2013-03-01
##
plot2020 <- function (summ, scenarios, designs, add = FALSE,
                      draw2020 = FALSE, draw1x = FALSE, ...) {
    if (missing(scenarios))
        scenarios <- 1:dim(summ$mean)[1]
    if (missing(designs))
        designs <- 1:dim(summ$mean)[3]
    recapt <- apply(summ$mean[scenarios,c('n','ncapt'),designs], c(1,3), diff)
    RSE.D <- summ$mean[scenarios, 'RSE.D', designs, drop=FALSE]

    if (!add)
        plot(recapt, RSE.D, xaxs = 'i', yaxs = 'i',
             xlim=c(0,1.1*max(recapt)), ylim=c(0,1.1*max(RSE.D)),
             xlab = 'Number of recaptures', ylab = 'RSE(D-hat)', type='n', ...)
    if (draw2020)
        lines(c(0,20,20),c(0.2,0.2,0), lty = 2)
    if (draw1x)
        lines(x <- 0:max(recapt), 1/x^0.5)
    points (recapt, RSE.D, ...)
    invisible(recapt)
}

## Examples
## plot2020(summary(outFM, sparse=T), scenarios = 1:8, pch=16)
## plot2020(summary(outFM, sparse=T), scenarios = 9:16, pch=1, add=T)
## mtext (side=3, line = 1, 'The 20:20 rule and 1/nrecapt^0.5', cex=0.8)
## legend (68,0.63, pch=c(1,16), legend=c("'male' g0=0.035, sigma=12.5 km",
##                                "'female' g0=0.013, sigma=6.5 km"), cex=0.9)
