
##==============================================================================
# shadowbox: Boxes with shadow
##==============================================================================

shadowbox <- function (box.type = "rect", mid, radx, rady=radx, shadow.size=0.01,
     shadow.col="grey", box.col="white", lcol="black", lwd=1, dr=0.01,
     angle=0, len=1, nr=5, rx=rady, theta = 90, ...)  {


  pin   <- par ("pin")                 # size of plotting region,  inches
  dd    <- c(shadow.size, -shadow.size*pin[1]/pin[2])     # scaled size of shadow

  if (box.type %in% c("rect", "square"))   {
    xy     <- cbind(c(mid[1]-radx, mid[1]-radx, mid[1]+radx, mid[1]+radx),
                    c(mid[2]-rady, mid[2]+rady, mid[2]+rady, mid[2]-rady))
    xyshad <- xy + matrix(nrow = 4, ncol = 2, byrow=TRUE, data=c(dd[1], dd[2]))
    if (angle != 0) {
      xy <- rotatexy  (xy, angle)
      xyshad <-rotatexy(xyshad, angle)
    }
    if (shadow.size>0)
      polygon(xyshad[, 1], xyshad[, 2], border=NA, col=shadow.col)
    polygon(xy[, 1], xy[, 2], lwd=lwd, col=box.col, border=lcol, ...)

  } else if (box.type %in% c("ellipse", "circle")) {
    if (shadow.size>0)
      filledellipse(mid=mid+dd, rx1=radx, ry1=rady, col=shadow.col,
                    dr=dr, angle=angle)
    filledellipse(mid=mid, rx1=radx, ry1=rady, col=box.col, dr=dr,
                  lwd=lwd, lcol=lcol, angle=angle, ...)

  } else if (box.type == "round") {
    if (shadow.size>0)
      roundrect(mid+dd, radx, rady, col=shadow.col, lcol=NA, dr=dr,
                rx=rx, angle=angle)
    roundrect(mid, radx, rady, lwd=lwd, col=box.col, lcol=lcol,
              rx=rx, angle=angle, ...)

  } else if (box.type =="diamond")            {
    xx    <- c(mid[1]-radx, mid[1], mid[1]+radx, mid[1])
    yy    <- c(mid[2], mid[2]+rady, mid[2], mid[2]-rady)
    xy    <- cbind(xx, yy)
    xyshad <- cbind(xx+dd[1], yy+dd[2])
    if (angle != 0) {
      xy <- rotatexy  (xy, angle)
      xyshad <-rotatexy(xyshad, angle)
    }
    if (shadow.size>0)
      polygon(xyshad[, 1], xyshad[, 2], border=NA, col=shadow.col)
    polygon(xy[, 1], xy[, 2], lwd=lwd, col=box.col, border=lcol, ...)

  } else if (box.type == "cylinder")         {
    if (shadow.size>0)
      filledcylinder(mid=mid+dd, rx=radx, ry=rady, len=len,
                     col=shadow.col, dr=dr, angle=angle)
    filledcylinder(mid=mid, rx=radx, ry=rady, len=len, col=box.col,
      dr=dr, lcol=lcol, lwd=lwd, angle=angle, ...)

  } else if (box.type== "hexa")             {
    if (shadow.size>0)
      filledmultigonal(mid=mid+dd, rx=radx, ry=rady, col=shadow.col,
        nr=6, angle=angle)
    filledmultigonal(mid=mid, rx=radx, ry=rady, col=box.col, lwd=lwd,
      nr=6, lcol=lcol, angle=angle, ...)

  } else if (box.type== "multi")             {
    if (shadow.size>0)
      filledmultigonal(mid=mid+dd, rx=radx, ry=rady, col=shadow.col,
        nr=nr, angle=angle)
    filledmultigonal(mid=mid, rx=radx, ry=rady, col=box.col, lwd=lwd,
      nr=nr, lcol=lcol, angle=angle, ...)

  } else if (box.type == "parallel") {
   # center box object
        shift <- 0.5*2*rady/tan(theta*pi/180)
        xy <- cbind(c(mid[1] - radx - shift,
                      mid[1] - radx - shift + 2*rady/tan(theta*pi/180),
                      mid[1] + radx - shift + 2*rady/tan(theta*pi/180),
                      mid[1] + radx - shift),
                      c(mid[2] - rady, mid[2] + rady,
                        mid[2] + rady, mid[2] - rady))
        xyshad <- xy + matrix(nrow = 4, ncol = 2, byrow = TRUE, data = c(dd[1], dd[2]))
        if (angle != 0) {
            xy <- rotatexy(xy, angle)
            xyshad <- rotatexy(xyshad, angle)
        }
        if (shadow.size > 0)
            polygon(xyshad[, 1], xyshad[, 2], border = NA, col = shadow.col)
        polygon(xy[, 1], xy[, 2], lwd = lwd, col = box.col, border = lcol,
            ...)
    }


  #else box.type=="none"
}
