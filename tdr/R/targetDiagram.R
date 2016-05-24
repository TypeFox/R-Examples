## `data` must include these columns: nrmse, nmbe, sdm, sdo

targetDiagram <- function(data, class = '',
                          xlab = expression("RMSEc"%.%"sign("*sigma^"*"*")"),
                          ylab = 'MBE',
                          auto.key = list(space = 'right'),
                          default.scales = list(cex = 0.6),
                          scales = list(),
                          type = 'quantile', cuts = seq(0.25, 1, .25),
                          par.settings = tdTheme(),
                          ...){
    
    data <- prepareData(data)
    circle <- makeCircles(data, type, cuts)
    radius <- unique(circle$r)

    ## Make formula
    ff <- as.formula(paste('nmbe ~ nrmsec * sign(difSD)',
                           ifelse(class == '', '', paste('|', class))
                           )
                     )
    ## Configure axis
    scales <- modifyList(default.scales, scales)
    ## Generate graphic
    xyplot(ff, 
           data = data,
           circle = circle,
           xlab = xlab, ylab = ylab,
           aspect='iso', scales = scales,
           xlim = extendrange(circle$x),
           ylim = extendrange(circle$y, f = 0.1),
           auto.key = auto.key,
           par.settings = par.settings,
           panel = function(..., circle){
               ## Vertical and Horizontal Axis
               panel.abline(h=0,v=0, lwd = 0.6, col='gray')
               ## Data
               panel.xyplot(...)
               ## Quantile circles
               panel.xyplot(circle$x, circle$y,
                            lwd = 0.6, type = 'l', col = 'darkgrey')
               ## Labels of circles
               panel.text(0, -radius, labels = signif(radius, 2),
                          pos = 1, offset = 0.05,
                          cex = scales$cex)
           }, ...)

}



