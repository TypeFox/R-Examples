setMethod("plot", signature(x = "abc"),
          function (x, xlim = c(0, ceiling(log(length(x@abc$Accum.Abund)))),
                    ylim = c(0, 100), yaxp = c(0, 100, 10), lty.bio = 'dotted',
                    lty.abu = 'solid', lwd = 2, col.bio = 'black', col.abu = 'black',
                    xlab = expression('Species Rank'~(Log[e]~Scale)),
                    ylab = 'Cumulative Dominance %', leg = TRUE, W = TRUE, col.pol = '#f5f5f5',
                    ...) {
            
            plot(x = log(1:length(x@abc$Accum.Abund)), y = x@abc$Accum.Abund, type = 'n',
                 xlim = xlim, ylim = ylim, ylab = ylab, xlab = xlab, axes = FALSE, ...)
            
            polygon(x = c(log(1:length(x@abc$Accum.Biomass)),
                          sort(log(1:length(x@abc$Accum.Abund)), decreasing=TRUE),
                          log(1:length(x@abc$Accum.Biomass))[1]),
                    y = c(sort(x@abc$Accum.Biomass, decreasing=FALSE),
                          sort(x@abc$Accum.Abund, decreasing=TRUE),
                          sort(x@abc$Accum.Biomass)[1]), col = col.pol, border = NA)
            
            lines(x = log(1:length(x@abc$Accum.Biomass)), y = sort(x@abc$Accum.Biomass),
                  lty = lty.bio, lwd = lwd, col = col.bio)
            
            lines(x = log(1:length(x@abc$Accum.Abund)), y = sort(x@abc$Accum.Abund,
                                                                 decreasing=FALSE), lty = lty.abu, lwd = lwd, col = col.abu)
            
            axis(side = 1, at = 0:ceiling(log(length(x@abc$Accum.Abund))),
                 labels = round(exp(0:ceiling(log(length(x@abc$Accum.Abund))))))
            
            axis(side = 2, yaxp = yaxp)
            
            if (W == TRUE)
              legend(x = 'topleft', legend = paste('W =', x@W.Stat[2], sep=' '),
                     bty = 'n')
            else (W == FALSE)
            legend(x = 'topleft', legend='', bty = 'n')
            
            if (leg == TRUE)
              legend(x = 'bottomright', legend = c('Biomass', 'Abundance'), lwd = lwd,
                     lty = c(lty.bio, lty.abu), bty = 'n', col = c(col.bio, col.abu))
            else (leg == FALSE)
            legend(x = 'bottomright', legend = '', bty = 'n')
          }
)