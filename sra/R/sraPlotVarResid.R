sraPlotVarResid <-
function (srafit, series = levels(srafit$data$rep)) 
{
    plot(NULL, NULL, axes = FALSE, xlim = c(0, max(srafit$data$gen)), 
        ylim = c(-max(abs(unlist(srafit$vresiduals))), max(abs(unlist(srafit$vresiduals)))), 
        xlab = "", ylab = "Residuals")
    axis(1, labels = FALSE, tick = TRUE)
    axis(2)
    abline(0, 0, lty = 3)
    for (r in names(srafit$pred)) {
        points(srafit$data$gen[srafit$data$rep == r], srafit$vresiduals[[as.character(r)]], 
            type = "b", lty = 2)
    }
}
