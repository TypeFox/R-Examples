ppccNormMultiplyCensored <-
function (x, censored, censoring.side, prob.method, plot.pos.con) 
{
    qqPlot.list <- qqPlotCensored(x = x, censored = censored, 
        censoring.side = censoring.side, prob.method = prob.method, 
        plot.pos.con = plot.pos.con, distribution = "norm", plot.it = FALSE)
    m.tilda <- qqPlot.list$x
    ss.m <- sum(m.tilda^2)
    c.vec <- m.tilda/sqrt(ss.m)
    index <- !qqPlot.list$Censored
    cor(c.vec[index], qqPlot.list$y[index])
}
