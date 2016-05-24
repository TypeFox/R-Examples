
##==============================================================================
# splitarrow: segmented arrow between several points
##==============================================================================

splitarrow <- function(from, to, lwd=2, lty=1, lcol="black", 
     arr.col=lcol, arr.side=2, arr.pos=0.5, centre=NULL, dd=0.5, ...)  {

  sarr <- function(p1, p2, drawarr)   {
    if (drawarr)
      m1<<-rbind(m1, straightarrow (from=p1, to=p2, arr.pos=arr.pos, lwd=lwd,
                     lty=lty, lcol=lcol, arr.col=arr.col, ...)) else
    segments(p1[1], p1[2], p2[1], p2[2], lwd=lwd, lty=lty, col=lcol)
  }

  m1   <- NULL
  From <- matrix(ncol=2, data=from)
  To   <- matrix(ncol=2, data=to  )

  meanFrom <- colMeans(From)
  meanTo   <- colMeans(To)

  if (is.null(centre))
    centre <- meanFrom+ dd*(meanTo-meanFrom)
  for (i in 1:nrow(From))
    sarr(From[i, ],  centre, 1 %in% arr.side)
  for (i in 1:nrow(To))
    sarr(centre, To[i, ]   , 2 %in% arr.side)

  splitarrow <- m1

}
