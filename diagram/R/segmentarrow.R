
##==============================================================================
# segmentarrorw: 3 segmented arrow between two points (left-vertical-right)
##==============================================================================

segmentarrow <- function(from, to, lwd=2, lty=1, lcol="black",
    arr.col=lcol, arr.side=2, arr.pos=0.5, path = "LVR", dd=0.5, ...)  {

  sarr <- function(p1, p2, drawarr) {
    if (drawarr)
      m1<<-rbind(m1, straightarrow (from=p1, to=p2, arr.pos=arr.pos, lwd=lwd,
                     lty=lty, lcol=lcol, arr.col=arr.col,...)) else
    segments(p1[1], p1[2], p2[1], p2[2], lwd=lwd, lty=lty, col=lcol)
  }
  m1 <- NULL
  if (is.null(path)) path <- ifelse (from[1]==to[1], "LVR", "UHD")
  if (path == "LVR") {dx <- -dd; dy <- 0}     # left vertical right
  if (path == "RVL") {dx <-  dd; dy <- 0}     # right vertical left
  if (path == "UHD") {dx <- 0  ; dy <- dd}    # up horizontal down
  if (path == "DHU") {dx <- 0  ; dy <- -dd}   # down horizontal up

  if (path %in% c("LVR", "RVL"))  {
    sarr(  from,                    c(from[1]+dx, from[2]   ), 1 %in% arr.side)
    sarr(c(from[1]+dx, from[2]   ), c(from[1]+dx,  to[2]+dy) , 2 %in% arr.side)
    sarr(c(from[1]+dx,  to[2]+dy) , to                       , 3 %in% arr.side)
  } else {
   sarr(  from,                   c(from[1]   , from[2]+dy), 1 %in% arr.side)
   sarr(c(from[1]   , from[2]+dy),c(to[1]     , from[2]+dy), 2 %in% arr.side)
   sarr(c(to[1]     , from[2]+dy),  to                     , 3 %in% arr.side)
  }
  segmentarrow <- m1       # coordinates (x,y) where arrowhead is drawn
}
