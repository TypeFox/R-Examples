`pairscor.fnc` <-
function(data, hist = TRUE, smooth = TRUE,
  cex.points = 1,  col.points = "darkgrey") {
  panel.hist <- function(x, ...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, ...)
  }

  pairscor.lower <- function(x, y, ...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
	  m = cor.test(x, y)
	  r = round(m$estimate, 2)
	  p = round(m$p.value, 4)
	  rtxt = paste("r =", r)
	  ptxt = paste("p =", p)

	  options(warn=-1)  # ignore warnings
	  m2 = cor.test(x, y, method="spearman")
	  r2 = round(m2$estimate, 2)
	  p2 = round(m2$p.value, 4)
	  rtxt2 = paste("rs =", r2)
	  ptxt2 = paste("p =", p2)
	  options(warn=0)

	  text(0.5, 0.8, rtxt)
	  text(0.5, 0.6, ptxt)
	  lines(c(0.2,0.8),c(0.5,0.5))
	  text(0.5, 0.4, rtxt2)
	  text(0.5, 0.2, ptxt2)
  }
  panel.smooth2 = function (x, y, col = par("col"), bg = NA, pch = par("pch"),
    cex = 1, span = 2/3, iter = 3, ...) {
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok))
        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),
            col = "black", ...)
  }


  if (hist == TRUE) {
    if (smooth == TRUE) {
	     pairs(data, 
         diag.panel = panel.hist,
         lower.panel = pairscor.lower, 
         upper.panel = panel.smooth2, col = col.points,
           cex = cex.points)
    } else {
       pairs(data, 
         diag.panel = panel.hist,
         lower.panel = pairscor.lower) 
    }
  } else {
    if (smooth == TRUE) {
	    pairs(data, lower.panel = pairscor.lower, 
        upper.panel = panel.smooth2, col = col.points,
        cex = cex.points)
    } else {
	    pairs(data, lower.panel = pairscor.lower) 
    }
  }
}

