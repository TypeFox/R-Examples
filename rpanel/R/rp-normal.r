rp.normal <- function(y, ylab = deparse(substitute(y)), 
                 panel.plot = TRUE, hscale = NA, vscale = hscale) {

   if (is.na(hscale)) {
      if (.Platform$OS.type == "unix") hscale <- 1
      else                             hscale <- 1.4
      }
   if (is.na(vscale)) 
      vscale <- hscale

   normal.draw <- function(panel) {
      panel$stan <- max(panel$stan, sd(panel$y) / 100)
      with(panel, {
         hist(y, prob = TRUE, col = "red", main = "", xlim = xlim, ylim = ylim)
         box()
         par.text <- character(0)
         xgrid <- seq(par()$usr[1], par()$usr[2], length = 100)
         if (curve.showing) {
            ygrid <- dnorm(xgrid, mean, stan)
            lines(xgrid, ygrid, col = "blue", lwd = 2)
            par.text <- paste("Mean =", as.character(signif(mean, 4)),
                              "  sd =", as.character(signif(stan,   4)))   
            title(par.text, col.main = "blue", line = 2, cex.main = 1)
            }
         if (fitted.showing) {
            ygrid <- dnorm(xgrid, mean(y), sd(y))
            lines(xgrid, ygrid, col = "green", lwd = 2)
            par.text <- paste("Mean =", as.character(signif(mean(y), 4)),
                              "  sd =", as.character(signif(sd(y),   4)))   
            title(par.text, col.main = "green", line = 1, cex.main = 1)
            }
         })
      panel
      }

   normal.redraw <- function(panel) {
     rp.tkrreplot(panel, plot)
     panel
     }

   hst <- hist(y, plot = FALSE)
   normal.panel  <- rp.control("Normal fitting", y = y, ylab = ylab,
                        ylim = c(0, max(hst$density) * 1.4),
                        xlim = c(min(y) - 2 * sd(y), max(y) + 2 * sd(y)),
                        mean = runif(1, 0.9 * min(y) + 0.1 * max(y), 0.7 * min(y) + 0.3 * max(y)),
                        stan = runif(1, 0.8 * sd(y), 1.4 * sd(y)),
                        curve.showing = FALSE, fitted.showing = FALSE)
   if (panel.plot & !require(tkrplot)) {
      warning("the tkrplot package is not available so panel.plot has been set to FALSE.")
      panel.plot <- FALSE
      }
   if (panel.plot) {
      rp.tkrplot(normal.panel, plot, normal.draw, pos = "right",
                 hscale = hscale, vscale = vscale, background = "white")
      action <- normal.redraw
      }
   else
      action <- normal.draw
   rp.checkbox(normal.panel, curve.showing, 
                        labels = "Show normal density", action = action)
   rp.doublebutton(normal.panel, mean, diff(range(y)) / 50, 
                        title = "Mean", action = action)
   rp.doublebutton(normal.panel, stan, sd(y) / 50,
                        title = "sd", action = action)
   rp.do(normal.panel, action)
   rp.checkbox(normal.panel, fitted.showing, 
                        labels = "Show fitted density", action = action)

   invisible()
   }
