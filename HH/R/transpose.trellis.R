transpose <- function(x)
  UseMethod("transpose")

transpose.trellis <- function(x) {
  if (length(dim(x)) != 2)
    return(rbind(x))
  y <- x
  ## dim is based on condlevels.  therefore new.order comes first before condlevels and index.cond
  new.order <- as.vector(t(array(1:prod(dim(x)), dim(x))))
  y$panel.args <- y$panel.args[new.order]
  y$condlevels <- rev(x$condlevels)
  y$index.cond <- rev(x$index.cond)
  y$layout <- rev(x$layout)
  y$par.settings$layout.heights$axis.panel <- x$par.settings$layout.widths$axis.panel
  y$par.settings$layout.widths$axis.panel <- x$par.settings$layout.heights$axis.panel
  y$par.settings$layout.heights$strip <- rev(x$par.settings$layout.widths$strip.left)
  y$par.settings$layout.widths$strip.left <- rev(x$par.settings$layout.heights$strip)
  if (is.list(y$x.limits)) {
    x.limits.all <- range(unlist(y$x.limits))
    for (ix in 1:length(y$x.limits))
      y$x.limits[[ix]] <- x.limits.all
  }
  if (is.list(y$y.limits)) {
    y.limits.all <- range(unlist(y$y.limits))
    for (iy in 1:length(y$y.limits))
      y$y.limits[[iy]] <- y.limits.all
  }
  y$x.between <- x$y.between
  y$y.between <- x$x.between
  if (is.list(y$x.scales))
    y$x.scales$at <- x$y.scales$at[new.order]
  if (is.list(y$y.scales))
    y$y.scales$at <- x$x.scales$at[new.order]
  y$packet.sizes <- t(x$packet.sizes)
  y
}

transpose.default <- function(x)
  t(x)

aperm.trellis <- function(a, perm, ...)  {
  y <- a
  if (is.null(a$layout))
    y$layout <- dim(y)
  if (length(y$layout) != length(dim(y)))
    stop("Please change layout to match dim.", call.=FALSE)
  ## dim is based on condlevels.  therefore new.order comes first before condlevels and index.cond
  new.order <- as.vector(aperm(array(1:prod(dim(y)), dim(y)), perm))
  y$panel.args <- y$panel.args[new.order]
  y$condlevels <- y$condlevels[perm]
  y$index.cond <- y$index.cond[perm]
  y$layout <- y$layout[perm]
  if (is.list(y$x.limits)) {
    x.limits.all <- range(unlist(y$x.limits))
    for (ix in 1:length(y$x.limits))
      y$x.limits[[ix]] <- x.limits.all
  }
  if (is.list(y$y.limits)) {
    y.limits.all <- range(unlist(y$y.limits))
    for (iy in 1:length(y$y.limits))
      y$y.limits[[iy]] <- y.limits.all
  }
  if (is.list(y$x.scales)) y$x.scales$at <- y$y.scales$at[new.order]
  if (is.list(y$y.scales)) y$y.scales$at <- y$x.scales$at[new.order]
  y$packet.sizes <- aperm(y$packet.sizes, perm)
  y
}

## cbind.trellis <- function(..., deparse.level=1,
##                           combineLimits=TRUE, useOuterStrips=TRUE) {
##   tmp <- rbind.trellis(..., deparse.level=deparse.level,
##                        combineLimits=FALSE, useOuterStrips=FALSE)
##   cdddots <- transpose(tmp)
##   if (useOuterStrips) cdddots <- useOuterStrips(cdddots)
##   if (combineLimits) cdddots <- combineLimits(cdddots)
##   cdddots
## }
