
plot.tclust <-
function (x,  ...)
{
  if (x$int$dim[2] == 1)
    .plot.tclust.1d (x, ...)
  else if (x$int$dim[2] == 2)
    .plot.tclust.2d (x, ...)
  else
    .plot.tclust.Nd (x, ...)
}

plot.tkmeans <-
function (x,  ...)
{
	plot.tclust (x, ...)
}

#######################
##  .plot.tclust.1d  ##
#######################

.plot.tclust.1d <-
function (x, xlab, ylab, xlim, ylim, tol = 0.95, tol.lwd = 1, tol.lty = 3, tol.col,
          jitter.y = FALSE, ...)
{
  if (x$int$dim[2] != 1)
    stop ("tclust object of dimension 1 expected.")
  
  if (is.null (x$par$x))
    stop ("dataset not included in tclust object - cannot plot object.")

  if (missing (xlab))
  {
    dn <- dimnames (x$par$x)
    if (is.null (dn[[2]]))
      xlab <- "x"
    else
      xlab <- dn[[2]][1]
  }

  if (missing (ylab) && jitter.y)
    ylab <- "(random jitter)"

  n <- x$int$dim[1]

  if (jitter.y)
    y <- runif (n, min = -1) / 4
  else
    y <- rep (0, n)

  if (is.numeric (tol) && length (tol) == 1 &&  0 < tol && tol < 1)
    tol.fact = sqrt (qchisq(tol, 1))
  else
    tol.fact <- NULL

  x.c = as.numeric (x$centers)
  if (!is.null (x$cov))		## tkmeans- objects don't have cov-info
  {
	  x.sd = sqrt (as.numeric (x$cov))

	  if (missing (xlim))
	  {
		xlim <- range (x$par$x)
		if (!is.null (tol.fact))
		  xlim <- range (xlim, x.c + x.sd * tol.fact, x.c - x.sd * tol.fact)
	  }
  }
  else
  {
	if (missing (xlim))
		xlim <- range (x$par$x)
  }

  if (missing (ylim))
	ylim <- c (-1, 1)

  X <- cbind (x$par$x, y)
  .plot.tclust.0 (x = x, X = X, xlab = xlab, ylab = ylab, axes = 1,
                 xlim = xlim, ylim = ylim, ...)
 
  .vline (x.c, 3, lty = 2, col = 1 + (1:x$k))    ##  cluster centers

	if (!is.null (tol.fact) &&
		!is.null (x$cov))		## tkmeans- objects don't have cov-info
  {
    #tol.fact = sqrt(qchisq(tol, 1))
    if (missing (tol.col))
      tol.col <- (1:x$k) + 1
    else
      tol.col <- rep (tol.col, x$k)
      
    tol.lty <- rep (tol.lty , x$k)
    tol.lwd <- rep (tol.lwd , x$k)   
    .vline (x.c + x.sd * tol.fact, 2, col = tol.col, lwd = tol.lwd,
            lty = tol.lty)
    .vline (x.c - x.sd * tol.fact, 2, col = tol.col, lwd = tol.lwd,
            lty = tol.lty)
  }
}

#######################
##  .plot.tclust.2d  ##
#######################

.plot.tclust.2d <-
function (x, xlab, ylab, tol = 0.95, tol.lwd = 1, tol.lty = 3, tol.col = 1,
          ...)
{
  if (nrow (x$centers) != 2)
    stop ("tclust object of dimension 2 expected.")

  if (is.null (x$par$x))
    stop ("dataset not included in tclust object - cannot plot object.")

  dn <- dimnames (x$par$x)
  if (is.list (dn) && length (dn[[2]]) == 2)
  {
    if (missing (xlab))
      xlab = dn[[2]][1]
    if (missing (ylab))
      ylab = dn[[2]][2]
  }
  else
  {
    if (missing (xlab))
      xlab = "x1"
    if (missing (ylab))
      ylab = "x2"
  }

  X <- cbind (x$par$x[, 1:2])
  .plot.tclust.0 (x = x, X = X, xlab = xlab, ylab = ylab, axes = 3, ...)

  if (!is.null (x$cov) && is.numeric (tol) && length (tol) == 1 &&  0 < tol && tol < 1)
  {
    tol.col <- rep (tol.col, x$k)
    tol.lty <- rep (tol.lty, x$k)
    tol.lwd <- rep (tol.lwd, x$k)

    tol.fact = sqrt(qchisq(tol, 2))  
    for (k in 1:x$k)
        .doEllipses (eigen = eigen (x$cov[,,k]), center = x$centers[,k],
        lwd = tol.lwd, lty = tol.lty[k], col = tol.col[k], size = tol.fact)
  }
}

#######################
##  .plot.tclust.Nd  ##
#######################

.plot.tclust.Nd <- 
function (x, xlab, ylab, ...)
{
  if (is.null (x$par$x))
    stop ("dataset not included in tclust object - cannot plot object.")

  if (missing (xlab))
    xlab <- "First discriminant coord."
  if (missing (ylab))
    ylab <- "Second discriminant coord."

  X <- discr_coords (x, x$par$equal.weights)

  .plot.tclust.0 (x = x, X = X, xlab = xlab, ylab = ylab, axes = 0, ...)
}

######################
##  .plot.tclust.0  ##
######################

.plot.tclust.0 <-
function (x, X, labels = c ("none", "cluster", "observation"), text,
          xlab, ylab, col, pch, by.cluster = TRUE, axes = 3, xlim, ylim, ...)
{

  if (by.cluster)
  {
    maxassig <- max (x$cluster)
    
    if (missing (col))
      col <- 1:(x$k+1)
    else
      col <- rep (col,  len = maxassig + 1 )  

    if (missing (pch))
      pch <- 1:(x$k+1)#rep (1, x$k+1)  #
    else
      pch <- rep (pch,  len = maxassig + 1 )
    col <- col[x$cluster + 1]    
    pch <- pch[x$cluster + 1]
  }
  else
  {
    if (missing (col))
      col <- x$cluster + 1
    if (missing (pch))
      pch <- 1
  }

  n <- x$int$dim[1]

  if (!missing (text))
    text <- rep (text, len = n)
  else if (!missing (labels))
  {
    labels <- match.arg(labels)
    if (labels == "cluster")
      text = paste (x$cluster)
    else if (labels == "observation")
      text = paste (1:nrow (X))
  }

  plot.new ()
  par (usr = .plot.tclust.calc.usr (X, xlim, ylim))

  if (missing (text))
    points (X[,1],X[,2], pch = pch, col = col)
  else
    text (X[,1], X[,2], labels = text, col = col)

  .plot.tclust.title (x, ...)

  axis.x <- axes %% 2        ## x axis 1 or 3
  axis.y <- axes >= 2        ## y axis 2 or 3

  cex <- par ("cex")
  if (!missing (xlab))
    mtext (side = 1, xlab, line = 1.5 + 1.5 * axis.x, cex = cex)
  if (!missing (ylab))
    mtext (side = 2, ylab, line = 1.5 + 1.5 * axis.y, cex = cex)

  box ()
  if (axis.x)
    axis (1)
  if (axis.y)
    axis (2)
}

##########################
##  .plot.tclust.title  ##
##########################

.plot.tclust.title <- function (x, main, main.pre, sub, sub1, ...)
{
  sub.par <- TRUE
  sub.restr <- FALSE

  sub.ovr <- missing (sub)
  if (!missing (sub) && is.character (sub))
  {
    sub.ovr <- TRUE
    if (sub == "/p")
      sub.par <- !(sub.restr <- FALSE)
    else if (sub ==  "/r")
      sub.par <- !(sub.restr <- TRUE)
    else if (sub ==  "/pr")
      sub.par <- sub.restr <- TRUE
    else
      sub.ovr <- FALSE
  }

	if (is.null (x$par$restr.C))
		sub.restr <- FALSE

  ralph <- round (x$par$alpha, 2)

  txt.par <- bquote(paste (k == .(x$par$k), ", ", alpha == .(ralph)))

  if (is.null(x$par$restr.C))
	txt.restr <- ""
  else if (x$par$restr.C == 2) ## this is the sigma - restriction, thus no restr.fact has to be printed.
    txt.restr <- paste ("restr = \"", x$par$restr, "\"", sep = "")
  else
    txt.restr <- paste ("restr = \"", x$par$restr, "\", restr.fact = ", x$par$restr.fact, sep = "")

  if (missing (main))
    main <- "Classification"  #"Cluster Assignment"
  else
    if (is.character (main))
      if (main == "/r" && !is.null (x$par$restr.C))
        main <- txt.restr
      if (main == "/p")
        main <- txt.par

  if (!missing (main.pre) && !is.null (main.pre))
    main <- paste (main.pre, main)

  if (sub.ovr)
    if (sub.par)
      sub  <- txt.par
    else if (!is.null (x$par$restr.C))
      sub <- txt.restr

  if (missing (sub1) && sub.restr && !(!sub.par && sub.ovr))
    sub1 <- txt.restr

  n.sub = .is.visible (sub) + .is.visible (sub1)    ##  number of subtitles to draw

  if (!is.null (main))
    title (main = main, line = ifelse (n.sub > 1, 2.3, 1.6))
  
  if (.is.visible (sub))
    mtext(sub, cex = 0.8, line = ifelse (n.sub > 1, 1.2, 0.3))

  if (.is.visible (sub1))
    mtext(sub1, cex = 0.8, line = ifelse (n.sub > 1, 0.1, 0.3))
}

.is.visible <- function (x)
{
  if (missing (x))
    return (FALSE)
  !is.null (x) && x != ""
}

.plot.tclust.calc.usr <- function (X, xlim, ylim, fact = 0.04)
{
  if (missing (xlim))
    xlim <- range (X[, 1])

  if (missing (ylim))
    ylim <- range (X[, 2])

  r <- cbind (xlim, ylim)
  rd <- apply (r, 2, diff)
#  r + rd * fact * c(-1,1)
  as.numeric (r + (c(-1, 1) %*% t (rd)) * fact)
}
