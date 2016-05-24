Showfry<-function(RDAT, shear= matrix(c(1, 1.2, 0,  1)), rad=75)
{

  if(missing(shear))
    {
      shr =0
      shear = matrix(c(1, 1.2, 0,  1), ncol=2)
    }
  if(missing(rad))
    {
      rad=75
    }
  
  x = RDAT$x
  y = RDAT$y
  APTS = cbind(x-mean(x),y-mean(y))
    
    NP  =  APTS %*% shear
    
    x = NP[,1]
    y = NP[,2]
    
    flag = sqrt( (x-mean(x))^2  + (y-mean(y))^2)<rad
    DAT = list(x=x[flag], y=y[flag])
    x = DAT$x
    y = DAT$y
    plot(x,y, asp=1, pch=21, col='red', bg='gold', ann=FALSE, axes=FALSE)
    
    FF = dofry(x, y )
    AF= plotfry(FF, dis=30)
    Z = xtractlip(AF)
    
    lines(cluster::predict.ellipsoid(Z$hull), col='red', lwd=1, lty=2)  
    points(Z$lip$x, Z$lip$y, col='green', pch=3, cex=.5)

  lines(Z$Elips,  col='brown', lwd=1.5)
 
  invisible(list(FF=FF, AF=AF, Z=Z))
  }
