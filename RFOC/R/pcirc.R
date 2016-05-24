`pcirc` <-
function(gcol='black', border='black', ndiv=36)
{
  if(missing(gcol)) { gcol='black' }
  if(missing(border)) { border='black' }
  if(missing(ndiv)) { ndiv=36 }


  phi = seq(0,2*pi, by=2*pi/ndiv);
  x = cos(phi);
  y = sin(phi);
  lines(x,y, col=border)
  lines(c(-1,1), c(0,0), col=gcol)
  lines(c(0,0), c(-1,1), col=gcol)

  

}

