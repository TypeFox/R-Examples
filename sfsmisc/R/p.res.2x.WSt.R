#### was part of ./p.goodies.R

### Exports :

### p.res.2x            Werner Stahels Plot; z.B Residuen gegen 2 x-Var.
### p.res.2fact         Aehnliche Idee: Residuen gegen 2 Faktoren (boxplots)

## p.wstPlot <- function(...)
## {
## warning("\n\n*** p.wstPlot(.) heisst neu p.res.2x(.)\n** Diese verwenden!\n")
##    p.res.2x(...)
## }

p.res.2x <- function(x, ...) UseMethod("p.res.2x")

p.res.2x.default <-
  function(x, y, z, restricted = NULL, size = 1, slwd = 1, scol = 2:3,
           xlab = NULL, ylab = NULL, main = NULL,
           xlim = range(x), ylim = range(y), ...)
{
  ## Purpose:  Stahels Residuen-Plot
  ## Author:   ARu , Date:  11/Jun/91
  ## Aenderungen: MMae, 30/Jan/92, Dez.94 --> help(p.res.2x)
  if(is.null(xlab)) xlab <- deparse(substitute(x))
  if(is.null(ylab)) ylab <- deparse(substitute(y))
  if(is.null(main)) main <- deparse(substitute(z))

  ok <- !(is.na(x) | is.na(y) | is.na(z))
  x <- x[ok]; y <- y[ok]; z <- z[ok]
  ##
  ##--- restrict z values: ---
  az <- abs(z)
  has.restr <-
    if(is.null(restricted)) FALSE else any(restr <- az > restricted)
  if(has.restr) {
    z[z >   restricted] <- restricted
    z[z < - restricted] <- - restricted
  }

  ##--- fix plot region: ---
  pcm <- par("pin") * 2.54              #damit in cm
  ##--- damit im Plot das Symbol wirklich die Groesse size hat:
  size <- size/(2 * sqrt(2))
  fx <- (size * diff(xlim))/(pcm[1] - 2 * size)/2
  fy <- (size * diff(ylim))/(pcm[2] - 2 * size)/2
  ##--
  plot(x, y, xlim = xlim + c(-1,1)* fx, ylim = ylim + c(-1,1)* fy, pch = ".",
       xlab = xlab, ylab = ylab, main = main, ...)

  ##--- draw symbols: ---
  z <- z/max(az, na.rm = TRUE)
  usr <- par("usr")
  sxz <-     diff(usr[1:2])/pcm[1] * size * z
  syz <- abs(diff(usr[3:4])/pcm[2] * size * z)
  if(length(scol) == 2) scol <- scol[1 + as.integer(z < 0)]
  segments(x - sxz, y - syz,  x + sxz, y + syz, lwd = slwd, col = scol)

  ##--- mark restricted observations: ---
  if(has.restr) {
    points((x - sxz)[restr], (y - syz)[restr], pch = 8, mkh = 1/40)
    points((x + sxz)[restr], (y + syz)[restr], pch = 8, mkh = 1/40)
  }
  invisible()
}

## graphics:::mosaicplot.formula as an example
p.res.2x.formula <- function(x = ~., data,
                             main = deparse(substitute(data)),
                             xlab = NULL, ylab = NULL, ...)
{
    ## Purpose:  plot residuals vs. two x's
    ## Author:   ARu , Date:  11/Jun/91
    ## Aenderungen: MMae, 30/Jan/92, Dez.94 / WSt
    ## --------------------------------------------------------------------------
    ## Arguments:
    ##   x         formula defining the variables zu be used, either
    ##             z ~ x + y
    ##             ~ x + y   in this case, data must inherit from  lm ,
    ##             and the residuals of  data  will be used as  z .
    ##   data      a data.frame or an  lm  or  aov  object.
    ##             In the latter case, p.res.2x will look for the data
    ##             that was used to fit the model.
    ##  restricted      absolute value which truncates the size.
    ##             The corresponding symbols are marked by stars.
    ##  size       the symbols are scaled so that 'size' is the size of
    ##             the largest symbol in cm.
    ##  slwd, scol line width and color to be used for the symbols
    ##  ...        additional arguments for the S-function 'plot'
    ## EXAMPLE :
    ## g.res2x(zz~.,data=data.frame(xx=rep(1:10,7),yy=rep(1:7, rep(10,7)),
    ##    zz=rnorm(70)), restr = 2, main = "i.i.d.  N(0,1) random residuals")
    ## --------------------------------------------------------------------------
    if(miss.main <- missing(main))
	force(main)
    formula <- as.formula(x)
    t.d <- if(inherits(data, "lm")) {
	if(miss.main) main <- paste0("residuals(", main, ")")
	if(!is.data.frame(t.d <- data$model)) {
	    ## try to look for the data that was used to fit the model.
	    cl <- data$call
	    i <- if("data" %in% names(cl)) "data" else 3 # try ..
	    t.d <- get(as.character(cl[[i]]))
	}
	if (length(formula) < 3) {
	    if(identical(format(formula), "~.")) formula <- formula(data)
	    formula <- update.formula(formula, residuals ~ .)
	    ## formula <- substitute(residuals ~ RHS, list(RHS = formula[[2]]))
	    cbind(t.d, residuals =  residuals(data))
	} else t.d
    } else
	data
    if (!is.data.frame(t.d)) {
	if(is.matrix(data)) data <- as.data.frame(data) else
	stop("data is not a data frame or 'lm' object with 'model' or existing data")
    }
    t.d <- na.omit(model.frame(formula, t.d))
    z <- t.d[,1]
    x <- t.d[,2]; if(is.null(xlab)) xlab <- names(t.d)[2]
    y <- t.d[,3]; if(is.null(ylab)) ylab <- names(t.d)[3]
    if(is.factor(x) && is.factor(y))
	p.res.2fact(x, y, z, main=main, xlab=xlab, ylab=ylab, ...)
    else {
	x <- as.numeric(t.d[,2])
	y <- as.numeric(t.d[,3])
	p.res.2x.default(x,y,z, main=main, xlab=xlab, ylab=ylab, ...)
    }
}

p.res.2fact <-
    function(x, y, z, restricted, notch = FALSE,
             xlab = NULL, ylab = NULL, main = NULL)
{
    if(is.null(xlab)) xlab <- deparse(substitute(x))
    if(is.null(ylab)) ylab <- deparse(substitute(y))
    if(is.null(main)) main <- deparse(substitute(z))

    ok <- !(is.na(x) | is.na(y) | is.na(z))
    x <- x[ok]; y <- y[ok]; z <- z[ok]
    x <- as.factor(x)
    y <- as.factor(y)
    lx <- levels(x);  ly <- levels(y)

    ##--- restrict z values: ---
    if(missing(restricted))  restr <- FALSE
    else {
        if(!is.numeric(restricted) || restricted <= 0)
            stop("'restricted' must be POSITIVE !")
        if(any(restr <- abs(z) > restricted)) {
            zorig <- z
            z[z >  restricted] <-   restricted
            z[z < -restricted] <- - restricted
        }
    }
    rz <- range(z)
    op <- par(mfrow = c(length(ly), 1), oma = c(5,6,6,0), mar = .1 + c(2,4,0,1))
    on.exit(par(op))
    for (yv in rev(ly)) {
        Ind <- y == yv
        plot (x[Ind], z[Ind], ylim = rz, xlab = "", ylab = yv, notch = notch)
        abline(h = 0, lty = 3, lwd = 0)
        if(any(II <- restr & Ind)) {
            ## boxplot creates a coord.system with x = [-4, 104]
            jx <- as.numeric(x[II])     #-- in 1:length(lx)..
            cat("..Cut z=",format(zorig[II])," at ",
                xlab,"=",x[II],",  ", ylab, "=",yv,"\n")
            points( u.boxplot.x(length(lx),jx) , z[II]*1.02, pch = 8, mkh = 1/25)
        }
    }
    mtext (xlab, side = 1, line = 1, outer = TRUE, cex = 1.3)
    mtext (ylab, side = 2, line = 3, outer = TRUE, cex = 1.3)
    mtext (main, side = 3, line = 2, cex = 1.5, outer = TRUE)
    if(any(restr)) message(sum(restr), " restricted observation(s)")
    invisible()
}


## Not sure if I want this (as global function).
## I had eliminated it long ago (from "SfS") but it's used above:

u.boxplot.x  <- function(n, j = 1:n, fullrange = 100)
{
  ## Purpose: Return the j-th x-coordinates in an 'n' side-by-side boxplot
  ## -------------------------------------------------------------------------
  ## Arguments: n : number of boxplots;  j: indices of boxplots
  ##  fullrange: x-coords as 'uniform' in [0,fullrange]  (f.=100, Splus 3.1,3.2)
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 19 Jan 95, 17:57
  cn <- fullrange/(3*n*(n+1))
  Dn <- cn*(3*n+2) ## Delta_{n}
  an <- cn*(2*n+1) ## a_{n}
  ## x(j) = an + (j-1)*Dn :
  an + (j-1)*Dn
}

