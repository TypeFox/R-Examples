
##==============================================================================
# treearrow: segmented arrow between several points
##==============================================================================

treearrow <- function(from, to, lwd=2, lty=1, lcol="black", 
    arr.col=lcol, arr.side=2, arr.pos=0.5, line.pos=0.5, path = "H", ...)  {

  sarr <- function(p1, p2, drawarr)  {
    if (drawarr)
      m1<<-rbind(m1, straightarrow (from=p1, to=p2, arr.pos=arr.pos, lwd=lwd,
                          lty=lty, lcol=lcol, arr.col=arr.col, ...)) else
    segments(p1[1], p1[2], p2[1], p2[2], lwd=lwd, lty=lty, col=lcol)
  }

  From <- matrix(ncol=2, data=from)
  To   <- matrix(ncol=2, data=to  )
  m1   <- NULL
  ifelse (path == "H",  ii <- 2, ii<-1)
  rF   <- range(From[, ii])
  rT <- range(To[, ii])

  ifelse (abs(min(rF)-max(rT)) <= abs(max(rF)-min(rT)),
     m2 <- min(rF) + line.pos * (max(rT) - min(rF)),
     m2 <- max(rF) + line.pos * (max(rF) - min(rT)) )


  if (path == "H")   {    # horizontal
    Line <- range(c(From[, 1], To[, 1]))
    segments(Line[1], m2, Line[2], m2, lwd=lwd, lty=lty, col=lcol)
    for (i in 1:nrow(From))
      sarr(From[i,],  c(From[i,1], m2), 1 %in% arr.side)
    for (i in 1:nrow(To)  )
      sarr(c(To  [i,1], m2), To[i,]   , 2 %in% arr.side)

  } else { # vertical
    Line <- range(c(From[,2], To[,2]))
    segments(m2, Line[1], m2, Line[2], lwd=lwd, lty=lty, col=lcol)
    for (i in 1:nrow(From))
      sarr(From[i,],  c(m2,From[i,2]), 1 %in% arr.side)
    for (i in 1:nrow(To)  )
      sarr(c(m2, To  [i,2]), To[i,]   , 2 %in% arr.side)
  }
  treearrow <- m1                  # coordinates (x,y) where arrowhead is drawn
}
