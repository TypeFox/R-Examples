# $Id: equivalence.xyplot.R,v 1.7 2005/10/10 10:14:43 andrewr Exp $

"equivalence.xyplot" <-
    function(formula,          
             alpha,      # The size; e.g. 0.05
             b0.ii,      # The half-length of the intercept r.s. (proportion)
             b1.ii,      # The half-length of the slope r.s. (proportion)
             add.smooth = FALSE,
             b0.absolute = FALSE,
             ...) {
      if (!b0.absolute) cat("Using b0.absolute = FALSE is deprecated and the argument will disappear in a future version.\n")
      bg <- trellis.par.get("background")
      bg$col <- "white"
      trellis.par.set("background", bg)
      print(xyplot(formula, subscripts = TRUE,
                   panel = function(x, y, subscripts, ...) {
                     x.bar <- mean(x, na.rm = TRUE)
                     the.model <- lm(y ~ x)
                     panel.xyplot(x, y, type = "n")
                     min.x <- min(x, na.rm = TRUE)
                     max.x <- max(x, na.rm = TRUE)
                     if (b0.absolute) y.poly <- x.bar + b0.ii * c(-1, 1, 1, -1)
                     else y.poly <- x.bar * (1 + b0.ii * c(-1, 1, 1, -1))
                     grid.polygon(x = c(min.x, min.x, max.x, max.x),
                                  y = y.poly,
                                  gp = gpar(col = "light gray",
                                    fill = gray(0.9)),
                                  default.units = "native",
                                  draw = TRUE)
                     panel.xyplot(x, y)
                     n <- sum(complete.cases(cbind(x, y)))
                     ybar.hat <- predict(the.model,
                                         newdata = data.frame(x = x.bar),
                                         se = TRUE)
                     t.quant <-  qt(1 - alpha/2, df.residual(the.model))
                     panel.arrows(x.bar, ybar.hat$fit +
                                  ybar.hat$se.fit * t.quant,
                                  x.bar, ybar.hat$fit -
                                  ybar.hat$se.fit * t.quant,
                                  col = "darkgrey", length = 0.05,
                                  angle = 90, code = 3)
                     se.slope <- coef(summary(the.model))[2, 2] 
                     panel.arrows(x.bar, ybar.hat$fit +
                                  se.slope * t.quant * x.bar,
                                  x.bar, ybar.hat$fit -
                                  se.slope * t.quant * x.bar, 
                                  col = "black", length = 0.05,
                                  angle = 90, code = 3)
                     panel.abline(a = coef(summary(the.model))[1, 1],
                                  b = 1 - b1.ii,
                                  col = "darkgrey", lty = 2)
                     panel.abline(a = coef(summary(the.model))[1, 1],
                                  b = 1 + b1.ii,
                                  col = "darkgrey", lty = 2)
                     panel.abline(coef(the.model))
                     if (add.smooth) panel.loess(x, y)
                   },
                   strip = strip.custom(style = 3), ...))
    }

