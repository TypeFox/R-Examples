#   Logistic regression cartoon

rp.logistic <- function(x, y, xlab = NA, ylab = NA, panel.plot = TRUE, panel = TRUE,
                hscale = NA, vscale = hscale, alpha = 0, beta = 0, 
                display = c("jitter" = FALSE, "regression line" = FALSE, 
                            "fitted model" = FALSE)) {

   if (is.na(hscale)) {
      if (.Platform$OS.type == "unix") hscale <- 1
      else                             hscale <- 1.4
      }
   if (is.na(vscale)) 
      vscale <- hscale
      
   if (is.na(xlab))
      xlab <- deparse(substitute(x))
   if (is.na(ylab))
      ylab <- deparse(substitute(y))
      
   if (is.matrix(y)) {
      n <- y[ , 1] + y[ , 2]
      y <- c(y[ , 1])
      }
   else {
      if (length(levels(factor(y))) > 2) 
         stop("y has more than two levels.")
      else {
         y <- as.numeric(factor(y)) - 1
         n <- rep(1, length(y))
         }
      }
      
   plot.logistic <- function(panel) {
   	
      if (!panel$panel.plot) {
         if (!panel$display["jitter"])
            yp <- panel$y
         else if (!panel$jitter.old)
            yp <- jitter(panel$y, amount = 0.02)
         else
            yp <- panel$y1
         panel$y1 <- yp
         panel$jitter.old <- panel$display["jitter"]
      }
      else
         yp <- panel$y1
      
      with(panel, {
         plot(x, yp / n, xlim = range(xgrid), xlab = xlab, ylab = ylab)
         abline(h = 0, lty = 3)
         abline(h = 1, lty = 3)
         title("Logistic regression", col.main = "red", line = 3, cex.main = 1)
         if (display["fitted model"]) {
            alpha1 <- par.est[1]
            beta1  <- par.est[2]
            clr    <- "green"
            lf     <- exp(alpha1 + beta1 * xgrid)
            lf     <- lf / (1 + lf)
            lines(xgrid, lf, lwd = 2, col = clr)
            if (alpha1 != 0) {
      	       int <- paste(signif(alpha1, 5))
      	       if (beta1 > 0) sgn <- " + " else sgn <- " - "
               }
            else {
               int <- ""
               if (beta1 > 0) sgn <- "" else sgn <- " - "
               }
            if (abs(beta1) != 1) slp <- paste(signif(abs(beta1), 5))
               else slp <- ""
            title(paste("log(p/(1-p)) = ", int, sgn, slp, " ", xlab, sep = ""), 
               col.main = clr, line = 1, cex.main = 1)            
            # title(paste(, signif(alpha1 - beta1 * mean(x), 5), "+", 
            #                            signif(beta1, 5), "x"), col = clr)
            }
         if (display["regression line"]) {
            alpha1 <- alpha - beta * mean(x)
            beta1  <- beta
            clr    <- "blue"
            lf     <- exp(alpha1 + beta1 * xgrid)
            lf     <- lf / (1 + lf)
            lines(xgrid, lf, lwd = 2, col = clr)
            if (alpha1 != 0) {
      	       int <- paste(signif(alpha1, 5))
      	       if (beta1 > 0) sgn <- " + " else sgn <- " - "
               }
            else {
               int <- ""
               if (beta1 >= 0) sgn <- "" else sgn <- " - "
               }
            if (abs(beta1) != 1) slp <- paste(signif(abs(beta1), 5))
               else slp <- ""
            title(paste("log(p/(1-p)) = ", int, sgn, slp, " ", xlab, sep = ""), 
               col.main = clr, line = 2, cex.main = 1)            
            # title(paste("log(p/(1-p)) =", signif(alpha1 - beta1 * mean(x), 5), "+", 
            #                            signif(beta1, 5), "x"), col = clr)
         }
      })
      panel
   }
      
   replot.logistic <- function(panel) {
      if (!panel$display["jitter"])
        panel$y1 <- y
      else if (!panel$jitter.old)
         panel$y1 <- jitter(panel$y, amount = 0.02)
      rp.control.put(panel$panelname, panel)
      rp.tkrreplot(panel, plot)
      panel$jitter.old <- panel$display["jitter"]
      panel
      }
            
   launch.lik <- function(panel) {
      with(panel, {
         data <- list(x = x, y = y, n = n)
         rp.likelihood(
              "sum((theta[1] + theta[2] * data$x) * data$y - log(1 + exp(theta[1] + theta[2] * data$x)) * data$n)", 
              data, par.est - 4 * par.se, par.est + 4 * par.se)
         })
      panel
      }

   dx      <- diff(range(x))
   xgrid   <- seq(min(x) - dx / 4, max(x) + dx / 4, length = 100)
   model   <- glm(cbind(y, n - y) ~ x, family = binomial)
   par.est <- coef(model)
   par.se  <- coef(summary(model))[,2]
   
   if (panel) {
      panel   <- rp.control(x = x, y = y, y1 = y, n = n, xgrid = xgrid, panel.plot = panel.plot,
                     alpha = alpha, beta = beta, par.est = par.est, par.se = par.se,
                     xlab = xlab, ylab = ylab, jitter.old = FALSE,
                     hscale = hscale, vscale = vscale, display = display)
      if (panel.plot & require(tkrplot)) {
         rp.tkrplot(panel, plot, plot.logistic, pos = "right", hscale = hscale, vscale = vscale,
                    background = "white")
         action <- replot.logistic
         }
      else
         action <- plot.logistic
      rp.doublebutton(panel, alpha, 0.05, action = action)
      rp.doublebutton(panel, beta,  3.5 / (diff(range(x)) * 50), action = action)
      rp.checkbox(panel, display, action, c("jitter", "regression line", "fitted model"))
      rp.button(panel, launch.lik, "Inspect log-likelihood")
      rp.do(panel, action)
   }
   else {
      panel <- list(x = x, y = y, y1 = y, n = n, xgrid = xgrid, panel.plot = FALSE,
                    alpha = alpha, beta = beta, par.est = par.est, par.se = par.se,
                    xlab = xlab, ylab = ylab, jitter.old = FALSE,
                    display = display)   
      plot.logistic(panel)
      invisible()
   }

   }
