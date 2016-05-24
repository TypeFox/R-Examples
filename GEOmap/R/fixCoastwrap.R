fixCoastwrap<-function(Z, maxdis=100)
{

  ##  this function assumes that the
######   polygons represent closed strokes.
###   the  points that wrap must come in pairs...
###   if not, skip
  
  if(missing(maxdis)) maxdis=100
  ww  = which( abs(abs(diff(Z$x))) > maxdis)

  nww = length(ww)

  if(nww<1) return(Z)

  
  
###   cbind(Z$x[ww], Z$x[ww+1])
###  print(paste("fixCoastwrap 1", nww))
### 
  
  if(RPMG::fmod(nww, 2)==1)
    {

   ###   print(paste("0 fixCoastwrap nww=", nww))
      u = par('usr')

      
      iw1 = which.min(Z$x-u[1])

      
      I2 = c( iw1:length(Z$x), 1:(iw1-1))
      
      Z = list(x=Z$x[I2], y=Z$y[I2])
      
      ww  = which( abs(abs(diff(Z$x))) > maxdis)
      
      nww = length(ww)
      
      if(nww>0)
        {
          ###  print(paste("1 fixCoastwrap nww=", nww))
          H = GEOmap.breakline(Z, ww)

          newx  = H$newx
          newy  = H$newy

        }
      else
        {

          newx  = Z$x
          newy  = Z$y
        ###  print(paste("2 fixCoastwrap nww=", nww))
          Z = list(x=c(min(Z$x) ,Z$x, max(Z$x))  ,  y=c(-89, Z$y, -89))
          return(Z)
        }


      
      if(nww==1)
        {
          Z = list(x=c(newx[[1]], NA, newx[[2]]), y=c(newy[[1]], NA, newy[[2]]) )
 
          # if(length(newx[[1]])>1 & length(newx[[2]])>1)
            {
          Z = list(x=c(min(Z$x) ,Z$x, max(Z$x))  ,  y=c(-90, Z$y, -90))
        }


          return(Z)
        }
      
      nww = length(ww)
    }


  j = 0
  newx = list()
  newy = list()

  if(nww<=1) return(Z)

  ## ww = c(0, ww, length(Z$x)+1)
####  build up a set of strokes separated by NA's
  
  H = GEOmap.breakline(Z, ww)

  newx=H$newx
  newy = H$newy

  
  if(FALSE)
    {
      lines(newx[[1]], newy[[1]], col='red')
      lines(newx[[2]], newy[[2]], col='green')
      lines(Z$x, Z$y, col='purple')


      polygon(c(u[1],Z$x, u[2])  , c(-90, Z$y, -90)  , col='purple')




      iw1 = which.min(Z$x-u[1])
      I = seq(from=1, to=length(Z$x), by=1)
      I2 = ((I-1)-iw1 )
      I2[I2<0] = I2[I2<0]

      u = par('usr')
      iw1 = which.min(Z$x-u[1])
      I2 = c( iw1:length(Z$x), 1:(iw1-1))
      Z = list(x=Z$x[I2], y=Z$y[I2])


    }

  if(nww==1)
    {


      return(list(x=c(newx[[1]], NA, newx[[2]]), y=c(newy[[1]], NA, newy[[2]]) ))




    }


####  if K is even there is a problem:
  K = length( newx)

  if(RPMG::fmod(K, 2)==0) {

   ### print(paste("even K fixCoastwrap 2", K) )

 ###    print(ww)
    return(list(x=c(newx[[1]], NA, newx[[2]]), y=c(newy[[1]], NA, newy[[2]]) ))
    

  }
  


  Side1.x = c(newx[[1]])
  
  Side1.y = c(newy[[1]])
  Side2.x = NULL
  Side2.y =  NULL

  ### print(paste("odd  K fixCoastwrap 2", K) )
  for(k in seq(from=2, to=K, by=2))
    {

      
      Side1.x = c(Side1.x, newx[[k+1]])
      Side1.y = c(Side1.y,  newy[[k+1]] )

      
      Side2.x = c(Side2.x, newx[[k]], NA)
      Side2.y = c(Side2.y,  newy[[k]] , NA)

      
    }




  ##  F.x = c(newx[j1:j2], newx[1:j3]  )  
  ##   F.y = c(newy[j1:j2], newy[1:j3]  )

  F.x =c(Side1.x, NA, Side2.x)
  F.y =c(Side1.y, NA, Side2.y)

  
  ##  polygon(F.x , F.y, col=rgb(.9,.9,1))

  return(list(x=F.x, y=F.y))

}

