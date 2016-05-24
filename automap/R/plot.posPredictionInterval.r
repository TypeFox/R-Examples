plot.posPredictionInterval = function(x, sp.layout = NULL, justPosition = TRUE, main = "Position prediction interval",...)
# Plot function for posPredictionInterval
{
    upper = automapPlot(x$pos_prediction,
                        zcol = "upper",
                        main = "Upper boundary",
                        sp.layout = sp.layout,
                        ...)

    lower = automapPlot(x$pos_prediction,
                        zcol = "lower",
                        main = "Lower boundary",
                        sp.layout = sp.layout,
                        ...)

    position2color = list(higher = "red", lower = "blue")
    position2color["not distinguishable"] = "yellow"
    color = sapply(levels(x$pos_prediction$position), FUN = function(x) position2color[[x]])
    position = spplot(x$pos_prediction,
                        zcol = "position",
                        main = main,
                        col.regions = color,
                        sp.layout = sp.layout,
						sub = paste(x$p, "% pred. interval relative to ", round(x$value, digits = 1), sep = ""),
                        ...)
    if(justPosition) print(position)
    else {                    
        print(lower, position = c(0,.5,.5,1), more = T)
        print(upper, position = c(.5,.5,1,1), more = T)
        print(position, position = c(0,0,.5,.5)) }
}
