sraPlotlegend <-
function (labels, estimates, AIC = NULL, confint = NULL, location = "topleft") 
{
    expr.labels <- sraFormatlegend(names = labels, values = estimates, 
        AIC = AIC, digits = 3)
    legend(x = location, as.expression(expr.labels), bg = "white")
}
