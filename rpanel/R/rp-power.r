rp.power <- function(panel = TRUE, panel.plot = TRUE, populations.showing = FALSE,
                     ngrid = seq(10, 300), mu1 = 0, mu2 = 1,
                     sigma = 1, n = 20, xgrid = seq(- 4, 5, length = 100),
                     popdens.lim = 0.7, hscale = 1, vscale = hscale) {

   powerplot.pars <- function(powerplot) {
      se              <- powerplot$sigma * 2 / sqrt(powerplot$ngrid)
      powerplot$pow   <- 1 - pnorm(2 - abs(powerplot$mu2 - powerplot$mu1) / se)
      rp.control.put(powerplot$panelname, powerplot)
      rp.do(powerplot, action.fn)
      powerplot
   }

   powerplot.draw <- function(powerplot) {
      powerplot$n <- max(powerplot$n, min(powerplot$ngrid))
      powerplot$n <- min(powerplot$n, max(powerplot$ngrid))
      with(powerplot, {
         if (populations.showing) par(mfrow = c(2, 1))
         par(mar = c(2, 3, 3.5, 1) + 0.1, mgp = c(1.0 , 0.2, 0), tcl = -0.2)
         plot(ngrid, pow, type = "n", ylim = c(0, 1), xlab = "Sample size, n", ylab = "Power")
         lines(ngrid, pow, col = "blue")
         lines(c(n, n, min(ngrid)), c(0, rep(pow[ngrid == n], 2)), lty = 2)
         title("Power tool", col.main = "red", line = 2.5)
         title(paste(  "mu1:",   format(round(mu1,   3), nsmall = 2), 
                     "  mu2:",   format(round(mu2,   3), nsmall = 2),
                     "  sigma:", format(round(sigma, 3), nsmall = 2)),
                     line = 1.5, cex.main = 0.85)
         title(paste("n:", n, "  Power:", format(round(pow[ngrid == n], 3), nsmall = 3)),
                     line = 0.5, cex.main = 0.85, col.main = "blue")
         if (populations.showing) {
            par(mar = c(2, 3, 1, 1) + 0.1, mgp = c(1.0, 0.2, 0), tcl = -0.2)
            plot(range(xgrid), c(0, popdens.lim), type = "n", 
                 xlab = "x", ylab = "Density")
            lines(xgrid, dnorm(xgrid, mu1, sigma))
            lines(xgrid, dnorm(xgrid, mu2, sigma))
            }
         par(mfrow = c(1, 1))
         })
      powerplot
   }

   powerplot.redraw <- function(panel) {
      rp.tkrreplot(panel, plot)
      panel
   }

   se  <- sigma * 2 / sqrt(ngrid)
   pow <- 1 - pnorm(2 - abs(mu2 - mu1) / se)
   
   if (panel) {
      if (panel.plot)
         action.fn <- powerplot.redraw
      else
         action.fn <- powerplot.draw
      power.panel <- rp.control("Power tool", populations.showing = populations.showing,
                        ngrid = ngrid, mu1 = mu1, mu2 = mu2,
                        sigma = sigma, n = n, xgrid = xgrid, pow = pow,
                        popdens.lim = popdens.lim, action.fn = action.fn)
      if (panel.plot)
         rp.tkrplot(power.panel, plot, powerplot.draw, pos = "right",
                    hscale = hscale, vscale = vscale, background = "white")
      rp.doublebutton(power.panel, n, 1, title = "n", action = action.fn)
      rp.doublebutton(power.panel, mu1, 0.01, title = "mu1", action = powerplot.pars)
      rp.doublebutton(power.panel, mu2, 0.01, title = "mu2", action = powerplot.pars)
      rp.doublebutton(power.panel, sigma, 0.01, title = "sigma", action = powerplot.pars)
      rp.checkbox(power.panel, populations.showing, labels = "Show populations",
         action = action.fn)
      rp.do(power.panel, powerplot.pars)
   }
   else {
   	   power.panel <- list(ngrid = ngrid, mu1 = mu1, mu2 = mu2,
                     sigma = sigma, n = n, xgrid = xgrid, pow = pow,
                     popdens.lim = popdens.lim)
   	   powerplot.draw(panel)
   }
   
   invisible()
}
