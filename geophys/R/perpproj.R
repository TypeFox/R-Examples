perpproj<-function( V1, V2, add=FALSE  )
  {
    if(missing(add)) add=FALSE
    ###  perpendicular projection routine 
###   get the projection points (altitude locations)
    ###  of the perpendicular legs of one point
    ###   extending down to another point on two vectors
####   
###  tails of each vector are assumed to be
    ### at point 1
    if(is.list(V1))
      {
        V1 = unlist(V1)

      }
    if(is.list(V2))
      {
        V2 = unlist(V2)

      }

#### V1 = c(      x1 ,      x2 ,      y1   ,    y2)
vlen = function(a){ return(sqrt(sum(a*a) )) }

###     need point of intersection of the two vectors
###    if they are not the same

    if(  !(V1[1]==V2[1] & V1[3] ==V2[3] ) )
      {

        ###  get intersecting point of 2 lines
      
       
            m1 = (V1[4]-V1[3])/(V1[2]-V1[1])
            m2 = (V2[4]-V2[3])/(V2[2]-V2[1])
            
            b1 = (V1[2]*V1[3]  -V1[1]* V1[4])/(V1[2]-V1[1])
            b2 =  (V2[2]*V2[3]  -V2[1]* V2[4])/(V2[2]-V2[1])
            
            x = -(b1-b2)/(m1-m2)
            y=m1*( x ) + b1

            V1 =c(x, V1[2], y, V1[4])
            
              V2 = c(x, V2[2], y, V2[4])
            
      }
    
    

v1 = c(  V1[2]-V1[1] ,    V1[4]-V1[3]      )
v2 = c(  V2[2]-V2[1] ,    V2[4]-V2[3]      )

    len1 = vlen(v1)
    len2 = vlen(v2)
    u1 = v1/len1
    u2 = v2/len2

  SENS1 = sign( AXB.prod(u1, u2)[3] )
 
    ####  if v1 or v2 are points and not vectors, need to do something different
VP = vecproj(v1, v2)

##   extend1 = sqrt(vlen(v2) ^2 - (VP$dis1)^2)

  POS1 = c(V1[1]+  VP$dis1 *u1[1],
  V1[3]+   VP$dis1*u1[2])
  ####    points(POS1[1], POS1[2], pch=6, col='blue')


  POS2 = c(V2[1]+  VP$dis2 *u2[1],
  V2[3]+   VP$dis2*u2[2])

    
  if(add)
    {
      points(POS1[1], POS1[2], pch=6, col='blue')
     points(POS2[1], POS2[2], pch=7, col='blue')
    }


    return(list(P1=POS1, P2=POS2))

  }
