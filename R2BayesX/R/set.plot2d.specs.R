set.plot2d.specs <-
function(nc, args, col.lines, is.bayesx)
{
  lwd <- args$lwd
  if(is.null(lwd) || any(is.na(lwd)))
    lwd <- rep(1, nc)
  else
    lwd <- rep(lwd, length.out = nc)
  lty <- args$lty
  if((is.null(lty) || any(is.na(lty))) && !is.bayesx[1L])
    lty <- rep(1, nc)
  if((is.null(lty) || any(is.na(lty))) && is.bayesx[1L])
    lty <- c(1, 0, 0, 0, 0)
  if(length(lty) == 1L || length(lty) < nc)
    lty <- rep(lty, length.out = nc)
  if(is.null(col.lines))
    col.lines <- rep("black", nc)
  else
    col.lines <- rep(col.lines, length.out = nc)
  args$lty <- lty
  args$lwd <- lwd
  args$col.lines <- col.lines

  return(args)
}

