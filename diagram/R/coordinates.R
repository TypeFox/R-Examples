
##==============================================================================
# coordinates: coordinates of components,  based on number of elements per row
##==============================================================================

coordinates <- function(pos=NULL, mx=0.0, my=0.0, N=length(pos),
      hor=TRUE, relsize=1)    {

  if (is.null(pos))  {  # positioned on a circle
    alpha0  <- pi/2
    alpha   <- alpha0 - (1:N) * 2 * pi/N
    xl      <- cos(alpha)*0.4 + 0.5
    yl      <- sin(alpha)*0.4 + 0.5
    elpos   <- cbind(xl, yl)

   } else {          # row- or column wise specification or fully specified
    nr    <- length(pos)

    if (length(pos) != N*2)    {   # row- or columnwise specification
      dx    <- 1/(max(pos))
      dy    <- 1/(nr)
      dr    <- min(dx, dy)

      elpos <-NULL
      ypos  <- rev(seq(dy*0.5, nr*dy, dy))

      for ( i in 1:nr) {
        dx   <- 1/(pos[i])
        xpos <- dx*0.5+ seq(0, pos[i]*dx, dx)
        for (j in 1:pos[i])
          elpos<-rbind(elpos, c(xpos[j], ypos[i]))
      }
      if (!hor)
        elpos <- rotatexy(elpos, angle=90, mid=c(0.5, 0.5))

    } else {                  # full specification
      elpos<-pos
    }
  }

  usr <- par("usr")
  if (relsize != 1) {
    dx    <- (usr[2]-usr[1])*(1-relsize)/2
    usr[1]<- usr[1]+dx
    usr[2]<-usr[2]-dx
    dy    <- (usr[4]-usr[3])*(1-relsize)/2
    usr[3]<-usr[3]+dy
    usr[4]<-usr[4]-dy

    }
  elpos[,1]  <-usr[1]+elpos[,1]*(usr[2]-usr[1])
  elpos[,2]  <-usr[3]+elpos[,2]*(usr[4]-usr[3])


  coordinates<-elpos + matrix(nrow = nrow(elpos), ncol = 2, byrow=TRUE, data=c(mx, my))
     # 2-columned matrix with coordinates (x, y) of each element
}

