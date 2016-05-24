HAMMERprojXY<-function(phi, lam)
   {
     ##  HAMMER Projection from Wikipedia
     ## phi = co-latitude (radians)
     ## lam = longitude (radians)

           sq2 = sqrt(2)
           cphi = cos(phi)
           sinlam = sin(lam/2)
           
           den = sqrt(1 + cphi * cos(lam/2))
           
           x = 2*sq2*cphi *sin(lam/2)/den
           y = sq2*sin(phi)/den

           return(list(x=x, y=y))

         }

