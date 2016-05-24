plotObsEst <-
function(yObs, yEst, delta, xlab = NULL,
ylab = NULL, title = NULL, legendplot = TRUE,
legendpos = "topleft", maxvalue = NULL, minvalue = NULL)
{
if(is.null(maxvalue))
maxvalue <- max(yObs, yEst)
if(is.null(minvalue))
minvalue <- min(yObs, yEst)
lim <- c(minvalue, maxvalue)
plot(yObs[delta == 1], yEst[delta == 1], pch = 1, xlim =
lim, ylim = lim, col = 2, main = title, xlab=xlab, ylab=ylab)
abline(0, 1)
points(yObs[delta == 0], yEst[delta == 0], pch = 4, col
= 4)
if(legendplot)
{
leg.txt <- c("Uncensor", "Censor")
legend(legendpos, leg.txt, pch = c(1,4), col = c(2, 4), cex = 3/4)
}
}
