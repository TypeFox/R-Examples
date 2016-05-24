DoMohr <-
function(Stensor=diag(c(3,2,1)) , axis=NULL )
{

  if(missing(axis)) { axis = NULL } 

  if(missing(Stensor))
    {
      Stensor=diag(c(3,2,1))

    }

####   s1  = Stensor

  di = dim(Stensor)

  if(di[1]==2)
    {
      ES = eigen(Stensor)

      sx = Stensor[1,1]
      sy = Stensor[2,2]
      txy = abs( Stensor[2,1] )


#### s1 = (sx+sy)/2 + sqrt(   ((sx-sy)/2)^2 + txy^2) 
#### s2  = (sx+sy)/2 - sqrt(   ((sx-sy)/2)^2 + txy^2) 

      taumax = mean(ES$values)
      ESave = mean(ES$values)

      thetap = atan2(2*txy, (sx-sy)) /2

      thetapd =thetap*180/pi

      thetas = atan2((sx-sy), -2*txy    )/2
      thetasd = thetas*180/pi



      Save = (sx+sy)/2


      sigman = Save+cos(2*thetap)*(sx-sy)/2+txy*sin(2*thetap)
      taunt = -(sx-sy)*sin(thetap)*cos(thetap) + txy*(cos(thetap)^2-sin(thetap)^2)
      
###Rmohr = sqrt( (sx-Save)^2  +  (txy)^2)


      ps1 = ES$values[1]
      ps2 = ES$values[2]

      Rmohr = ps1-Save

      ex = Save
      why = 0

      cmohr = GEOmap::darc( rad=Rmohr, ang1=0, ang2=360, x1=ex, y1=why, n=1)

      RNGM = range( cmohr$x)
      PXrange = c(0 , RNGM[2]+0.05*diff(RNGM), RNGM[1]-0.05*diff(RNGM))


      plot(range( PXrange ) , range(cmohr$y), type='n', asp=1, axes=FALSE, ann=FALSE, xpd=TRUE)

      grid(col=grey(.85))
      abline(v=0, h=0, lty=1, lwd=1.5)

      
      
      lines(cmohr$x, cmohr$y, lwd=2, col='brown' )

      points(ex, why)

      axis(1)
      axis(2)
      abline(v=ex, lty=3, col=grey(.75) )


      points(sx, -txy)

      points(sy, txy)

      points(ps1, 0, col='red')
      text(ps1, 0, labels=expression(sigma[1]), adj=c(1,1)   )

      text(sx,-txy, labels=expression(sigma[x]), adj=c(0,1)   )
      text(sy, txy, labels=expression(sigma[y]), adj=c(1,-0.5)   )



      points(ps2, 0, col='red')
      text(ps2, 0, labels=expression(sigma[2]), adj=c(1,1)   )

      segments(sx, -txy ,  sy, txy)

     

      mohrleg(ES)
  u = par('usr')
       
    tex1 = substitute(theta[p]==x , list(x=thetapd) )
        mtext(tex1, side = 1, line = -1, at=u[2], adj=1 )

    tex1 = substitute(rho==x , list(x=Rmohr) )
        mtext(tex1, side = 1, line = -2, at=u[2], adj=1 )



      
      out = list(sx=sx, sy=sy, txy=txy, s1=ps1, s2=ps2, ES=ES)

    }
  if(di[1]==3)
    {

      ######   get eigen value decomposition
      ES = eigen(Stensor)

      ####  extract eigen values for further work
      s1 = ES$values[1]
      s2 = ES$values[2]
      s3 = ES$values[3]

      ####  get loci of centers of circles
      x1 = mean(c(s1,s3))
      x2 = mean(c(s1,s2))
      x3 = mean(c(s2,s3))

      ####  get radii of three circles
      Rmohr1 = s1-x1
      Rmohr2 = s1 - x2
      Rmohr3 = abs(s2-x3)

      ###   calculate the points of the circles
      cmohr1 = GEOmap::darc( rad=Rmohr1, ang1=0, ang2=360, x1=x1, y1=0, n=1)
      cmohr2 = GEOmap::darc( rad=Rmohr2, ang1=0, ang2=360, x1=x2, y1=0, n=1)
      cmohr3 = GEOmap::darc( rad=Rmohr3, ang1=0, ang2=360, x1=x3, y1=0, n=1)

      #####  for plotting get the range of the plot
      RNGM = range( c(cmohr1$x,cmohr2$x, cmohr3$x    ))
      PXrange = c(0 , RNGM[2]+0.05*diff(RNGM), RNGM[1]-0.05*diff(RNGM))

      RNGY = range(c(cmohr1$y,cmohr2$y, cmohr3$y    ))

      ####  initialize the plot
      plot(range( PXrange) , RNGY, type='n', asp=1, axes=FALSE, ann=FALSE, xpd=TRUE)

      if(!is.null(axis))
        {
          for(iax in 1:length(axis)) { axis(axis[iax]) }


        }
      
      ####   fill circles, big circle is shaded, smaller ones are white
      polygon(cmohr1$x, cmohr1$y, col=grey(.9))
      polygon(cmohr2$x, cmohr2$y, col="white")
      polygon(cmohr3$x, cmohr3$y, col="white")

      ### plot grid
      
      grid(col=grey(.75))

      ##   add (x,y) or (sigma and tau)  axes
      abline(v=0, h=0, lty=1, lwd=1.5)

    ###   plot principle stresses
      points(c(s1, s2, s3), rep(0, 3), col='red')

      ##  draw Mohr's circles
      lines(cmohr1$x, cmohr1$y, lwd=2, col='brown' )
      lines(cmohr2$x, cmohr2$y, lwd=2, col='brown' )
      lines(cmohr3$x, cmohr3$y, lwd=2, col='brown' )

      ####   add the text indicatring the 3 principle stresses
      text(s1, 0, labels=expression(sigma[1]), adj=c(1,1)   )
      text(s2, 0, labels=expression(sigma[2]), adj=c(1,1)   )
      text(s3, 0, labels=expression(sigma[3]), adj=c(1,1)   )

      ####   add a little legend on upper right corner
      mohrleg(ES)
      
      ###  return stress tensor and eigen values/eigen vectors
     # out = ES
      out = list(Sin=Stensor,  ES=ES)
    }



  return(out)

}

