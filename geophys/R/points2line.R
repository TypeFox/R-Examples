  points2line<-function(Lp, VL )
      {
####  Lp is the location of a point, x,y
####  basepoint is the first point of VL
####  VL is a set of lines as a matrix


        vy = cbind(Lp$x-VL[ ,1] , Lp$y-VL[ ,2])

        v = cbind(VL[ ,3] - VL[ ,1] , VL[ ,4] - VL[ ,2] )
        dv = sqrt(v[, 1]^2+v[, 2]^2)
        dvy = sqrt(vy[, 1]^2+vy[, 2]^2)

        cv = (v[, 1]*vy[, 1] + v[, 2]*vy[,2])/(dv*dvy)

        rat =  dvy*cv/dv
        srat = sqrt( dvy^2 -  (dvy*cv)^2 )
        return(list(rat=rat, srat=srat))
      }
    
