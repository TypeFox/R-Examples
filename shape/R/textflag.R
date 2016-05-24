
textflag <- function(mid, radx, rady, rx = rady, dr = 0.01, 
                     col = femmecol(100), lcol = "white", 
                     bcol = lcol,    # color to remove ellipse
                     lwd = 2, angle = 0, lab = NULL, leftright = TRUE, 
                     tcol = NULL, ...) {
  wx <- rady * 2
  wy <- radx * 2
  if (is.null(rx)) rx <- rady
 
  # rectangle
 if (leftright)
  filledrectangle(wx = wx, wy = wy, col = col, 
                   mid = mid, angle = angle -90)
 else
  filledrectangle(wx = wy, wy = wx, col = col, 
                   mid = mid, angle = angle)

  # ellipse left to be cut away
  leftell <- getellipse(ry = rady, rx = rx, mid = c(mid[1] - 
         radx+rx, mid[2]), dr = dr, from = pi/2, to = 3/2 * pi)

  lx <- range(leftell[, 1])
  ly <- range(leftell[, 2])
  LL <- rbind(leftell, c(lx[1], ly[1]), c(lx[1], ly[2]), c(lx[2], ly[2]))       
  if (angle != 0)
    LL <- rotatexy(LL, angle = angle, mid = mid)

  polygon(LL, col = bcol, border = bcol)         

  # ellipse left to be cut away
  rightell <- getellipse(ry = rady, rx = rx, mid = c(mid[1] + 
         radx-rx, mid[2]), dr = dr, from = -pi/2, to = pi/2)

  lx <- range(rightell[, 1])
  ly <- range(rightell[, 2])
  RR <- rbind(rightell, c(lx[2], ly[2]), c(lx[2], ly[1]), c(lx[1], ly[1]))    
  if (angle != 0)
    RR <- rotatexy(RR, angle = angle, mid = mid)
    
  polygon(RR, col = bcol, border = bcol)         

  # lines  
  Lines <- rbind(leftell,rightell,leftell[1,])
  if (angle != 0)
    Lines <- rotatexy(Lines, angle = angle, mid = mid)
  
  lines(Lines, col = lcol, lwd = lwd)
  if (! is.null(lab)) text (mid[1], mid[2], lab, srt = angle, col = tcol, 
    ...)
}

 