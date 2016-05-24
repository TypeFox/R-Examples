normcheck = function(x, ...){
  UseMethod("normcheck")
}

normcheck.default = function(x, xlab = NULL, main = xlab, col = NULL, shapiro.wilk = FALSE, ...){
  
  if(is.null(xlab))
    xlab  =   deparse(substitute(x))

  if(is.null(col))
      col = "light blue"

  ## only grab the parameters that are going to be set
  oldPar = par(c("mfrow", "xaxs", "yaxs", "pty"))
  
  

  ## change the layout of the plotting window only if it has not already been set
  if(all(oldPar$mfrow == c(1, 1))){
    par(mfrow = c(1, 2), xaxs = "i", yaxs = "i", pty = "s")
    on.exit(par(mfrow = oldPar))
  }

  mx = mean(x)
  sx = sd(x)
  qqp = qqnorm(x)
  abline(c(mx,sx))
  
  
  if(shapiro.wilk){
    stest = shapiro.test(x)
    txt = paste("Shapiro-Wilk normality test", "\n",
               "W = ", round(stest$statistic,4), "\n",
               "P-value = ", round(stest$p.value,3), sep="")
    text(sort(qqp$x)[2], 0.99*sort(qqp$y)[length(qqp$y)], txt, adj = c(0,1))
  }
  
  h = hist(x, plot = FALSE)
  rx = range(x)
  xmin = min(rx[1],mx-3.5*sx,h$breaks[1])
  xmax = max(rx[2],mx+3.5*sx,h$breaks[length(h$breaks)])


  ymax = max(h$density,dnorm(mx,mx,sx))*1.05

  hist(x, prob = TRUE, ylim = c(0, ymax), xlim = c(xmin, xmax),
       xlab = xlab, col = col,
       main = main)
  box()

  x1 = seq(xmin,xmax, length = 100)
  y1 = dnorm(x1, mx, sx)
  lines(x1,y1, lwd = 1.5, lty = 3)
}

normcheck.lm = function(x, xlab = NULL, main = xlab, col = NULL, shapiro.wilk = FALSE, ... ){	
  if (missing(x) || (class(x) != "lm") )
		stop ("missing or incorrect lm object")

  
  if(is.null(xlab))
    xlab = paste("Residuals from lm(", as.character(x$call[2]),")", sep="")
  
  x = residuals(x)
	normcheck (x, xlab = xlab, main = main, col = col, shapiro.wilk = shapiro.wilk, ...)
}

