`PLTcirc` <-
function(gcol='black', border='black', ndiv=36, angs=c(-pi, pi), PLOT=TRUE, add=FALSE)
{
  if(missing(gcol)) { gcol='black' }
  if(missing(border)) { border='black' }
  if(missing(ndiv)) { ndiv=36 }
  if(missing(angs)) {  angs=c(-pi, pi)  }
  if(missing(PLOT)) {  PLOT=TRUE  }
  if(missing(add)) {  add=TRUE  }


  phi = seq(angs[1], angs[2], by=(angs[2]-angs[1])/ndiv);
  x = cos(phi);
  y = sin(phi);
  if(add==FALSE) { plot(x,y, type='n', ann=FALSE, axes=FALSE, asp=TRUE) }
  if(PLOT==TRUE)
    {
      lines(x,y, col=border)
      lines(c(-1,1), c(0,0), col=gcol)
      lines(c(0,0), c(-1,1), col=gcol)
    }
  invisible(list(x=x, y=y, phi=phi))

}

