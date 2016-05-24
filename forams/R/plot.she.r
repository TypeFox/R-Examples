setMethod("plot", signature(x = "she"), 
          function (x, pch = 20, pcol = 'black', pcex = 1, pbg = 'black', lcol = 'black',
                    lwd = 1, lty = 'dotted', ylab = expression('ln'~italic(E)), bty = 'l', ...) {
            
            if (is.vector(x@bi$N) == TRUE)
              vari = expression('ln'~italic(M))
            else
              vari = expression('ln'~italic(L))
            
            plot(log(x@bi[[4]]), log(x@bi$E), type = 'n', bty = bty, xlab = vari,
                 ylab = ylab, ...)
            
            lines(log(x@bi[[4]]), log(x@bi$E), lty = lty, lwd = lwd, col = lcol)
            
            points(log(x@bi[[4]]), log(x@bi$E), pch = pch, col = pcol, bg = pbg, cex = pcex)
          }
)