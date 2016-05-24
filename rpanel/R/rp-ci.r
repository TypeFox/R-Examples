#   Simulations of confidence intervals

rp.ci <- function(mu = 0, sigma = 1, sample.sizes = c(30, 50, 100, 200, 500), confidence = 0.95,
                  panel = TRUE, panel.plot = TRUE, hscale = NA, vscale = hscale) {

   if (is.na(hscale)) {
      if (.Platform$OS.type == "unix") hscale <- 1
      else                             hscale <- 1.4
      }
   if (is.na(vscale)) 
      vscale <- hscale

   ci.sim <- function(panel) {
      n     <- as.numeric(panel$ssize)
      mu    <- as.numeric(panel$pars[1])
      sigma <- as.numeric(panel$pars[2])
      X     <- matrix(rnorm(n * 100, mu, sigma), ncol = n)
      Xmean <- apply(X, 1, mean)
      Xsd   <- apply(X, 1, sd)
      limit <- qt(1 - (1 - as.numeric(panel$confidence)) / 2, n - 1) * Xsd / sqrt(n)
      lower <- Xmean - limit
      upper <- Xmean + limit
      panel$cover    <- ((lower < mu) & (mu < upper))
      panel$coverage <- panel$coverage + length(which(panel$cover))
      panel$nsim     <- panel$nsim + 100
      panel$lower    <- lower
      panel$upper    <- upper
      if (panel$nopanel) 
         ci.plot(panel)
      else      
         rp.control.put(panel$panelname, panel) 
      if (!panel$first) {
         if (panel$panel.plot) rp.tkrreplot(panel, plot)
            else               rp.do(panel, ci.plot)
         }
      panel$first <- FALSE
      panel
      }
      
   ci.plot <- function(panel) {
      n     <- as.numeric(panel$ssize)
      mu    <- as.numeric(panel$pars[1])
      sigma <- as.numeric(panel$pars[2])
      plot(mu + c(-5, 5) * sigma / sqrt(30), c(1, 100), type = "n", xlab = "y", ylab = "")
      colour         <- rep("blue", 100)
      colour[!panel$cover] <- "red"
      segments(panel$lower, 1:100, panel$upper, 1:100, col = colour)
      abline(v = mu, lty = 2)
      title("Simulated confidence intervals", col.main = "red", line = 3, cex.main = 1)
      title(paste("Observed coverage: ", round(100 * panel$coverage / panel$nsim, 2), "%", sep = ""), 
              col.main = "blue", line = 2, cex.main = 1)
      title(paste("(Simulation size: ", panel$nsim, ")", sep = ""), line = 1, cex.main = 1)
      panel
      }
      
   reset.coverage <- function(panel) {
      panel$coverage <- 0
      panel$nsim     <- 0
      rp.control.put(panel$panelname, panel)
      rp.do(panel, ci.sim)
      panel
      }

   if (panel.plot & !require(tkrplot)) {
      warning("the tkrplot package is not available so panel.plot has been set to FALSE.")
      panel.plot <- FALSE
      }
      
   if (panel) {
      panel <- rp.control("Simulated confidence intervals", panel.plot = panel.plot,
                    pars = c("mean" = mu, "s.d." = sigma), ssize = 30, confidence = confidence,
                    coverage = 0, nsim = 0, first = TRUE, nopanel = !panel)
      rp.do(panel, ci.sim)
      if (panel.plot)
         rp.tkrplot(panel, plot, ci.plot, pos = "right", hscale = hscale, vscale = vscale,
                    background = "white")
      rp.textentry(panel, pars, labels = c("mean", "s.d."), action = reset.coverage)
      rp.radiogroup(panel, ssize, sample.sizes, title = "Sample size", action = reset.coverage)
      rp.radiogroup(panel, confidence, c(0.90, 0.95, 0.99), title = "Confidence",
                    action = reset.coverage)
      rp.button(panel, title = "Sample", action = ci.sim, repeatdelay = 500, repeatinterval = 200)
      rp.button(panel, title = "Reset coverage count", action = reset.coverage)
      if (!panel.plot)
         rp.do(panel, reset.coverage)
   }
   else {
   	  panel <- list(pars = c("mean" = mu, "s.d." = sigma), ssize = 30, confidence = confidence,
                    coverage = 0, nsim = 0, first = TRUE, nopanel = !panel)
      ci.sim(panel)
   }
   
   invisible()
   }
