cd_plot <- function(x, ...) {
  UseMethod("cd_plot")
}

cd_plot.formula <- function(formula, data = list(),
  plot = TRUE, ylab_tol = 0.05,
  bw = "nrd0", n = 512, from = NULL, to = NULL,
  main = "", xlab = NULL, ylab = NULL, margins = c(5.1, 4.1, 4.1, 3.1),
  gp = gpar(), name = "cd_plot", newpage = TRUE, pop = TRUE,
  return_grob = FALSE, ...)
{
  ## extract x, y from formula
  mf <- model.frame(formula, data = data)
  if(NCOL(mf) != 2) stop("`formula' should specify exactly two variables")
  y <- mf[,1]
  if(!is.factor(y)) stop("dependent variable should be a factor")
  x <- mf[,2]
  if(!is.numeric(x)) stop("explanatory variable should be numeric")

  ## graphical parameters
  if(is.null(xlab)) xlab <- names(mf)[2]
  if(is.null(ylab)) ylab <- names(mf)[1]

  ## call default interface
  cd_plot(x, y,
    plot = plot, ylab_tol = ylab_tol,
    bw = bw, n = n, from = from, to = to,
    main = main, xlab = xlab, ylab = ylab, margins = margins,
    gp = gp, name = name, newpage = newpage, pop = pop, ...)
}

cd_plot.default <- function(x, y,
  plot = TRUE, ylab_tol = 0.05,
  bw = "nrd0", n = 512, from = NULL, to = NULL,
  main = "", xlab = NULL, ylab = NULL, margins = c(5.1, 4.1, 4.1, 3.1),
  gp = gpar(), name = "cd_plot", newpage = TRUE, pop = TRUE,
  return_grob = FALSE, ...)
{
  ## check x and y
  if(!is.numeric(x)) stop("explanatory variable should be numeric")
  if(!is.factor(y)) stop("dependent variable should be a factor")
  ny <- length(levels(y))

  ## graphical parameters
  if(is.null(xlab)) xlab <- deparse(substitute(x))
  if(is.null(ylab)) ylab <- deparse(substitute(y))
  if(is.null(gp$fill)) gp$fill <- gray.colors(ny)
  gp$fill <- rep(gp$fill, length.out = ny)

  ## unconditional density of x
  dx <- if(is.null(from) & is.null(to)) density(x, bw = bw, n = n, ...)
          else density(x, bw = bw, from = from, to = to, n = n, ...)
  x1 <- dx$x

  ## setup conditional values
  yprop <- cumsum(prop.table(table(y)))
  y1 <- matrix(rep(0, n*(ny-1)), nrow = (ny-1))

  ## setup return value
  rval <- list()

  for(i in 1:(ny-1)) {
    dxi <- density(x[y %in% levels(y)[1:i]], bw = dx$bw, n = n, from = min(dx$x), to = max(dx$x), ...)
    y1[i,] <- dxi$y/dx$y * yprop[i]
    rval[[i]] <- approxfun(x1, y1[i,], rule = 2)
  }
  names(rval) <- levels(y)[1:(ny-1)]

  ## use known ranges
  y1 <- rbind(0, y1, 1)
  y1 <- y1[,which(x1 >= min(x) & x1 <= max(x))]
  x1 <- x1[x1 >= min(x) & x1 <= max(x)]

  ## plot polygons
  if(plot) {
    ## setup
    if(newpage) grid.newpage()
    pushViewport(plotViewport(xscale = range(x1), yscale = c(0, 1),
      default.units = "native", name = name, margins = margins, ...))

    ## polygons
    for(i in 1:(NROW(y1)-1)) {
     gpi <- gp
     gpi$fill <- gp$fill[i]
     grid.polygon(x = c(x1, rev(x1)), y = c(y1[i+1,], rev(y1[i,])), default.units = "native", gp = gpi)
    }

    ## axes
    grid.rect(gp = gpar(fill = "transparent"))
    grid.xaxis()
    grid.yaxis(main = FALSE)
    equidist <- any(diff(y1[,1]) < ylab_tol)
    yat <- if(equidist) seq(1/(2*ny), 1-1/(2*ny), by = 1/ny) else (y1[-1,1] + y1[-NROW(y1), 1])/2
    grid.text(x = unit(-1.5, "lines"), y = unit(yat, "native"), label = levels(y),
      rot = 90, check.overlap = TRUE)

    ## annotation
    grid.text(xlab, y = unit(-3.5, "lines"))
    grid.text(ylab, x = unit(-3, "lines"), rot = 90)
    grid.text(main, y = unit(1, "npc") + unit(2, "lines"), gp = gpar(fontface = "bold"))

    ## pop
    if(pop) popViewport()
  }

  ## return conditional density functions
  if (plot && return_grob)
      invisible(structure(rval, grob = grid.grab()))
  else
      invisible(rval)

}
