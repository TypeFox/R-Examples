##
## dsplot.R - desirability plot
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.tu-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.tu-dortmund.de>
##

dsplot <- function(expr, f,
                   from=NULL, to=NULL, n=101,
                   show.zero=TRUE, interest=NULL,
                   main="Desirability Plot", sub=NULL,...) {
  sexpr <- substitute(expr)
  if (is.name(sexpr)) {
    fcall <- paste(sexpr, "(x)")
    expr <- parse(text = fcall)
  } else {
    if (!(is.call(sexpr) && match("x", all.vars(sexpr), nomatch = 0))) 
      stop("'expr' must be a function or an expression containing 'x'")
    expr <- sexpr
  }

  ## Evaluate expression:
  x <- seq(from, to, length.out=n)
  y <- eval(expr, envir=list(x=x), enclos=parent.frame())

  if (length(interest) > 0) {
    xi <- interest
    yi <- eval(expr, envir=list(x=xi), enclos=parent.frame())
    di <- f(yi)
  }

  ## Save par and restore on exit
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))

  ## Allocate space for main/sub title
  oma <- c(0, 0, 0, 0)
  if (!is.null(main)) oma[3] <- 2
  if (!is.null(sub)) oma[1] <- 2
  par(oma=oma)

  ## Layout of plots: 1/3 desirability, 2/3 expression
  layout(matrix(c(1, 2), ncol=2),
         widths=c(2/3, 1/3),
         heights=c(1,1))

  ## Expression plot
  par(mar=c(5, 2, 1, 3))
  plot(x, y, type="l", yaxt="n", ...)
  axis(4)
  ## Global y axis label
  mtext(sexpr, side=2, line=1)
  
  if (show.zero)
    abline(h=0, col="grey", lty=2)
  if (length(interest) > 0) {
    xmax <- par("usr")[2]
    ymin <- par("usr")[3]
    segments(xi, ymin, xi, yi, col="red", lty=2)
    segments(xmax, yi, xi, yi, col="red", lty=2)
  }
  ## Desirability plot
  par(mar=c(5, 0, 1, 2))
  y <- y[is.finite(y)] ## Remove +/- Inf values (poles)
  y <- seq(min(y), max(y), length.out=n)  
  plot(f(y), y,
       xlab=substitute(d(sexpr)), ylab=sexpr,
       type="l", yaxt="n", ...)
  axis(2, labels=FALSE)

  if (show.zero)
    abline(h=0, col="grey", lty=2)
  if (length(interest) > 0) {
    xmin <- par("usr")[1]
    segments(xmin, yi, di, yi, col="red", lty=2)
  }
  
  ## Add main/sub title to plot:
  if (!is.null(main))
    mtext(main, side=3, outer=TRUE)
  if (!is.null(sub))
    mtext(sub, side=1, outer=TRUE)
  invisible()
}
