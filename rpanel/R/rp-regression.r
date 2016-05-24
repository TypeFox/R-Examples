#   Simple regression with R graphics

rp.regression <- function (x, y, ylab = NA, x1lab = NA, x2lab = NA, xlab = NA, 
           panel = TRUE, panel.plot = TRUE, hscale = NA, vscale = hscale,
           model = "None", line.showing = TRUE, residuals.showing = FALSE, size = 3, col = "red") {

   if (is.na(hscale)) {
      if (.Platform$OS.type == "unix") hscale <- 1
      else                             hscale <- 1.4
      }
   if (is.na(vscale)) 
      vscale <- hscale

rp.regression1 <- function(x, y, ylab, xlab, panel.plot, hscale = NA, vscale = hscale) {
                         	
   scatter.draw <- function(object) {
      with(object, {
         clr <- rep(1, length(x))
         if (!any(is.na(coords))) {
            if (Influence["move points horizontally"]) x[ind] <- coords[1]
            if (Influence["move points vertically"])   y[ind] <- coords[2]
            if (Influence["move points horizontally"] | Influence["move points vertically"])
               clr[ind] <- 2
         }
         plot(x, y, xlim = range(xlim, x), ylim = range(ylim, y), 
                 col = clr, xlab = xlab, ylab = ylab)
         title("Simple linear regression", col.main = "red", line = 3, cex.main = 1)
         line.text <- character(0)
         if (Display["regression line"]) {
            intercept.adj <- intercept - slope * (mean(x))
            abline(intercept.adj, slope, col = "blue")
            line.text <- paste(line.text,
                            "Intercept =", as.character(signif(intercept.adj, 4)),
                            "    Slope =", as.character(signif(slope,         4)))
            if (Display["residuals"]) {
               fitted.values <- intercept.adj + slope * x
               if (!Display["fitted line"]) segments(x, y, x, fitted.values, col = "red")
               ss <- syy - sxy^2 / sxx + length(y) * (mean(y) - intercept)^2 +
                      sxx * (slope - sxy / sxx)^2
               ss.text <- paste("        SS =", signif(ss, 5))
            }
            else
               ss.text <- ""
            if (intercept.adj != 0) {
         	   int <- paste(signif(intercept.adj, 5))
         	   if (slope > 0) sgn <- " + " else sgn <- " - "
            }
            else {
         	   int <- ""
         	   if (slope > 0) sgn <- "" else sgn <- " - "
            }
            if (abs(slope) != 1) slp <- paste(signif(abs(slope), 5))
            else                 slp <- ""
            title(paste("E(", ylab, ") = ", int, sgn, slp, " ", xlab, ss.text, sep = ""), 
               col.main = "blue", line = 2, cex.main = 1)
         }
         if (Display["fitted line"]) {
            model <- lm(y ~ x)
            cfs   <- coef(model)
            abline(cfs[1], cfs[2], col = "green")
            if (Display["residuals"]) {
               segments(x, y, x, fitted(model), col = "red")
               ss.text <- paste("       RSS =", signif(sum((y - fitted(model))^2), 5))
            }
            else
               ss.text <- ""
            if (cfs[1] != 0) {
      	       int <- paste(signif(cfs[1], 5))
      	       if (cfs[2] > 0) sgn <- " + " else sgn <- " - "
            }
            else {
               int <- ""
               if (cfs[2] > 0) sgn <- "" else sgn <- " - "
            }
            if (abs(cfs[2]) != 1) slp <- paste(signif(abs(cfs[2]), 5))
               else slp <- ""
            title(paste("E(", ylab, ") = ", int, sgn, slp, " ", xlab, ss.text, sep = ""), 
               col.main = "green", line = 1, cex.main = 1)            
            # title(line.text, col = "red")
         }
      })
      object$parplt <- par("plt")
      object$parusr <- par("usr")
      object
      }

   find.pt <- function(panel, x, y) {
	  tol.x <- diff(range(panel$x)) / 25
	  tol.y <- diff(range(panel$y)) / 25
	  d.pts <- (panel$x - x)^2 + (panel$y - y)^2
	  panel$ind <- min(which(d.pts == min(d.pts)))
	  if ((abs(panel$x[panel$ind] - x) < tol.x) & 
	      (abs(panel$y[panel$ind] - y) < tol.y))
         panel$coords <- c(x, y)
      else
         panel$coords <- rep(NA, 2)
      rp.control.put(panel$panelname, panel)
      rp.tkrreplot(panel, plot)
      panel
   }
   
   drag <- function(panel, x, y) {
	  if (!any(is.na(panel$coords)))
         panel$coords <- c(x, y)
      rp.control.put(panel$panelname, panel)
      rp.tkrreplot(panel, plot)
      panel
   }
   
   release <- function(panel, x, y) {
      panel$coords  <-  rep(NA, 2)
      rp.control.put(panel$panelname, panel)
      rp.tkrreplot(panel, plot)
      panel
   }
   
   replot.regression1 <- function(panel) {
      rp.tkrreplot(panel, plot)
      panel
   }
   
   intercept.initial <- runif(1, 0.7*min(y) + 0.3*max(y), 0.3*min(y) + 0.7*max(y)) 
   intercept.delta   <- diff(range(y)) / 200
   slope.delta       <- (diff(range(y)) / diff(range(x))) / 50
   # slope.delta       <- abs(coef(lm(y ~ x))[2]) / 50
   sxx               <- sum((x - mean(x))^2)
   syy               <- sum((y - mean(y))^2)
   sxy               <- sum((x - mean(x)) * (y - mean(y)))
   xlim              <- range(x)
   ylim              <- range(y)

   scatter.panel <- rp.control("Scatterplot",
                      x = x, y = y, xlab = xlab, ylab = ylab,
                      sxx = sxx, syy = syy, sxy = sxy, xlim = xlim, ylim = ylim,
                      intercept = intercept.initial, slope = 0,
                      Display = c("regression line" = line.showing,
                         "residuals" = residuals.showing, "fitted line" = FALSE),
                      coords = rep(NA, 2))
   if (panel.plot) {
      rp.tkrplot(scatter.panel, plot, scatter.draw, find.pt, drag, release, 
               hscale = hscale, vscale = vscale, pos = "right", background = "white")
      action.fn <- replot.regression1
      }
   else
      action.fn <- scatter.draw
   rp.doublebutton(scatter.panel, intercept, intercept.delta,
                      repeatinterval = 20,
                      title = "Intercept", action = action.fn)
   rp.doublebutton(scatter.panel, slope, slope.delta,
                      repeatinterval = 20,
                      title = "Slope     ", action = action.fn)
   rp.checkbox(scatter.panel, Display, action.fn, 
                      c("regression line", "residuals", "fitted line"))
   if (panel.plot)
      rp.checkbox(scatter.panel, Influence, action.fn, 
                      c("move points horizontally", "move points vertically"))
   rp.do(scatter.panel, action.fn) 
   invisible(scatter.panel)
}

   if (is.na(ylab)) ylab <- deparse(substitute(y))
   x.name <- deparse(substitute(x))
   if (is.vector(x)) {
      if (is.na(xlab)) xlab <- x.name
      rp.regression1(x, y, ylab, xlab, panel.plot, hscale = hscale, vscale = vscale)
   }
   else if (is.matrix(x)) {
   	  x.names <- dimnames(x)[[2]]
      name.comp<-if (!is.null(x.names) & !all(x.names == "")) x.names
                 else {if (!is.null(attributes(x)$names)) attributes(x)$names
                       else outer(x.name, c("[1]", "[2]"), paste, sep = "")}
      if (is.na(x1lab)) x1lab <- name.comp[1]
      if (is.na(x2lab)) x2lab <- name.comp[2]
   	  rp.regression2(y, x[ , 1], x[ , 2], ylab, x1lab, x2lab, panel, model, 
   	                    residuals.showing, size, col)
   }
      
}
