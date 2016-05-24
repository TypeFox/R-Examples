pointproj<-function(P1, VEC )
  {
    ####   calculates the perpendicular distance of point P1 to line VEC
  if(is.list(P1))
      {
        P1 = cbind(P1$x, P1$y)

      }
    if(is.list(VEC))
      {
        VEC = c(VEC$x, VEC$y)

      }

 
  N = length(P1[,1])
  
    ### if P1 is a vector matrix or points get the 

   ##  vector pointing from v1 to v2
    v = c(VEC[2]-VEC[1], VEC[4]-VEC[3])
    Lv = sqrt(v[1]^2+v[2]^2)

  ####  the length of the line is zero, return null
    if(Lv<=0) { return(NULL) }
  
  #### to origin of first point in vector
  vex =   cbind( P1[ ,1] - VEC[1], P1[ ,2] - VEC[3])
  ##  arrows(VEC[1], VEC[3], P1[ ,1], P1[ ,2])

###  set up cos and sin of rotation
  cosTHETA = (VEC[2]-VEC[1])/Lv
  sinTHETA = (VEC[4]-VEC[3])/Lv

###  set up rotation matrix
  M = matrix( c(cosTHETA, sinTHETA, -sinTHETA, cosTHETA), ncol=2, byrow=TRUE)
###  rotate vectors
 H2 =  M %*% t(vex)
#####   zero out y-part
 H3 = cbind(H2[1, ] , rep(0, N))
####   unrotate  
 H4 =   solve(M) %*% t(H3)
###   translate back to correct position
EX = H4[1,]+ VEC[1]
WHY = H4[2,]+VEC[3]

  dis = abs(H2[2, ])
  
    return(list(x=EX, y=WHY, d=dis ))
  }
