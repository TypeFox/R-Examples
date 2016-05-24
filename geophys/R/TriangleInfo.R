TriangleInfo<-function(P1, P2=c(0,1), P3=c(1,0) , add=FALSE)
  {

    
    ###  P1  P2 and P3 are 3 points in a plane (2D)
    
   ####   Planer calculations
   ####  given a triangle (three non-colinear points, return info that is useful

    if(missing(P2) | missing(P3))
      {
        ###########  if P2 and P3 is missing then P1 is a matrix (3 by 2 ) with the three points

        if(is.matrix(P1))
          {
            P2 = P1[2,]
            P3 = P1[3,]
            P1 = P1[1,]
          }
        if(is.list(P1))
          {
            P2 = c(P1$x[2], P1$y[2])
            P3 = c(P1$x[3], P1$y[3])
            P1 = c(P1$x[1], P1$y[1])

          }

        
      }

    
######  check for colinearity 
    X = rbind(P1, P2, P3)
    X = cbind(X, rep(1, 3))
    DX = det(X)
    if(DX==0)
      {
        cat("POINTS ARE COLINEAR!!")
        return(NULL)
      }
###################################
    

    if(missing(add)) { add=FALSE }

    
    if(is.list(P1)) P1 = unlist(P1)
    if(is.list(P2)) P2 = unlist(P2)
    if(is.list(P3)) P3 = unlist(P3)


    medcol = 'green'
    altcol = 'cyan'
    Abiscol = 'purple'
    PerpBiscol = 'blue'
#################   add=TRUE
    

####   set vector from each point to the next
    v3 = c( P2[1]-P1[1], P2[2]-P1[2] )
    v1 = c( P3[1]-P2[1], P3[2]-P2[2] )
    v2 = c( P1[1]-P3[1], P1[2]-P3[2] )
    
####   get sense of triangle (clock wise or counter clockwise)
    SENS1 = sign( AXB.prod(v3, -v2)[3] )

###  centroid

    
    CEN = c( (P1[1]+P2[1]+P3[1])/3 , (P1[2]+P2[2]+P3[2])/3 )

   if(add) points(CEN[1], CEN[2], pch=3)


    
    if(add){
      arrows(P1[1], P1[2], P1[1]+v3[1] , P1[2]+v3[2] , length=.1 )
       arrows(P2[1], P2[2], P2[1]+v1[1] , P2[2]+v1[2]  , length=.1 )
      arrows(P3[1], P3[2], P3[1]+v2[1] , P3[2]+v2[2]  , length=.1 )
    }

    
########  these are the lengths of each side.
    ####   the lengths are opposite the designatedP1 corner
    
    a =  sqrt( (v1[1])^2 + (v1[2])^2 )
    b  = sqrt( (v2[1])^2 + (v2[2])^2 )
    c  = sqrt( (v3[1])^2 + (v3[2])^2 )


    cang2 = sum( v1*(-v3) )/(a*c)
    angrad2 = acos(cang2)
    angdeg2 = angrad2*180/pi

    cang3 = sum((-v1)*v2)/(a*b)
    angrad3 = acos(cang3)
    angdeg3 = angrad3*180/pi


    cang1= sum((v3)*(-v2))/(b*c)
    angrad1 = acos(cang1)
    angdeg1 = angrad1*180/pi


    P = a+b+c
    s = P/2

    K = sqrt(s*(s-a)*(s-b)*(s-c))

########    K = (a*c*sin(angrad1))/2
########    K = (b*c*sin(angrad3))/2
########  sqrt(s*(s-a)*(s-b)*(s-c))

    ###   h's are the altitudes = the perpendicular distances to each leg
    ####   from the opposite point
    ha = 2*K/a
    hb = 2*K/b
    hc = 2*K/c

    ########   the m's are the lengths of the median from points to opposite sides
    ma = sqrt(2*b^2+2*c^2-a^2)/2
    mb = sqrt(2*a^2+2*c^2-b^2)/2
    mc = sqrt(2*a^2+2*b^2-c^2)/2

    #####   the t's are the lengths of the bisectors 
    ta = sqrt(b*c*(1 - (a^2/(b+c)^2)))
    tb = sqrt(a*c*(1 - (b^2/(a+c)^2)))
    tc = sqrt(a*b*(1 - (c^2/(a+b)^2)))


    z = (-v1[2])  * v2[1] - v2[2] * (-v1[1])

    Area = abs(z)
    
    r = K/s

    R  = (a*b*c)/(4*K)
    


########## directions of the medians

    U1 = c( (P2[1]+P3[1])/2 ,   (P2[2]+P3[2])/2)
    UM1 = c( U1[1]-P1[1] , U1[2]-P1[2])
   if(add)  arrows(P1[1], P1[2], U1[1], U1[2], col=medcol , length=.1)

    U2 = c( (P1[1]+P3[1])/2 ,   (P1[2]+P3[2])/2)
    UM1 = c( U2[1]-P2[1] , U2[2]-P2[2])
  if(add)   arrows(P2[1], P2[2], U2[1], U2[2], col=medcol , length=.1)

    U3 = c( (P1[1]+P2[1])/2 ,   (P1[2]+P2[2])/2)
    UM1 = c( U3[1]-P3[1] , U3[2]-P3[2])
  if(add)   arrows(P3[1], P3[2], U3[1], U3[2], col=medcol , length=.1)

    MEDs = rbind(U1, U2, U3)

########## directions of the altitudes



    perpc = c(-1*SENS1*v3[2],SENS1*v3[1])
    perpb = c(-1*SENS1*v2[2],SENS1*v2[1])
    perpa = c(-1*SENS1*v1[2],SENS1*v1[1])

    ####  these are the direction vectors
    perpc = perpc/sqrt(sum(perpc^2))
    perpb = perpb/sqrt(sum(perpb^2))
    perpa = perpa/sqrt(sum(perpa^2))
    
    ####  these are the points
    Hc = c(P3[1]-hc*perpc[1], P3[2] -hc*  perpc[2])
    Hb = c(P2[1]-hb*perpb[1] , P2[2] -hb*  perpb[2])
    Ha = c(P1[1]-ha*perpa[1], P1[2] -ha*  perpa[2])


    ALTs = rbind(Ha,Hb,Hc)


    
    
    
   if(add)  arrows(P3[1], P3[2],  Hc[1], Hc[2] , col=altcol  , length=.1)
    
   if(add)  arrows(P2[1], P2[2], Hb[1], Hb[2]  , col=altcol  , length=.1)
    
   if(add)  arrows(P1[1], P1[2], Ha[1], Ha[2]  , col=altcol  , length=.1)
    



    
    ####  intersection of the perpendicular bisectors, OH


    

####  find intersection point

    ab = sqrt( (Hb[1]-P3[1])^2+(Hb[2]-P3[2])^2 )
    cd = hc

    HRAT =  ab*b/cd 

    IH = c(P3[1]+HRAT* (Hc[1]-P3[1])/hc  ,  P3[2]+HRAT* (Hc[2]-P3[2])/hc)

  if(add)   points(IH[1], IH[2], col=altcol)
  if(add)   points(CEN[1], CEN[2], col=altcol)
    
################    directions of the angular Bisectors
    B1 = SENS1*angrad1/2
    BI1 = c(v3[1]*cos(B1)-v3[2]*sin(B1) , v3[1]*sin(B1)+v3[2]*cos(B1))
    BI1 = BI1/sqrt(sum(BI1*BI1))
###   DI1 = c
    DI1 = sin(angrad2)*c/sin((pi)-angrad2-angrad1/2)
    
    B1 = SENS1*angrad2/2
    BI2 = c(v1[1]*cos(B1)-v1[2]*sin(B1) , v1[1]*sin(B1)+v1[2]*cos(B1))
    BI2 = BI2/sqrt(sum(BI2*BI2))
    DI2 = sin(angrad3)*a/sin((pi)-angrad3-angrad2/2)
    
 ####   if(add) arrows( P2[1],P2[2], P2[1]+DI2*BI2[1]  , P2[2]+DI2*BI2[2],  col=Abiscol  , length=.1  ) 

    B1 = SENS1*angrad3/2
    BI3 = c(v2[1]*cos(B1)-v2[2]*sin(B1) , v2[1]*sin(B1)+v2[2]*cos(B1))
    BI3 = BI3/sqrt(sum(BI3*BI3))
    DI3 = sin(angrad1)*b/sin((pi)-angrad1-angrad3/2)
 ####   if(add) arrows( P3[1],P3[2], P3[1]+DI3*BI3[1]  , P3[2]+DI3*BI3[2],  col=Abiscol   , length=.1 )


    
    AngBis = cbind(
      c(P1[1]+DI1*BI1[1], P2[1]+DI2*BI2[1] , P3[1]+DI3*BI3[1]),
      c(P1[2]+DI1*BI1[2], P2[2]+DI2*BI2[2],   P3[2]+DI3*BI3[2]))
    

    if(add)
      {
        arrows( P1[1],P1[2],AngBis[1,1]   ,AngBis[1,2] ,  col=Abiscol   , length=.1)
        arrows( P2[1],P2[2], AngBis[2,1]   ,AngBis[2,2]   ,  col=Abiscol  , length=.1  ) 
        arrows( P3[1],P3[2], AngBis[3,1]   ,AngBis[3,2],  col=Abiscol   , length=.1 )

      }
    

    V2 = list(x=c(P2[1],P2[1]+a*BI2[1] ), y=c(P2[2], P2[2]+a*BI2[2] )) 
    V1 =  list(x=c(P3[1], P3[1]+b*BI3[1]  ), y=c(P3[2], P3[2]+b*BI3[2]) )

        ###  get intersection of 2 vectors
    CIRCinside = unlist(Sect2vex(V1, V2))
 
################
   if(add)
     {
       dcirc =  GEOmap::darc(rad = r, ang1 = 0, ang2 = 360, x1 = CIRCinside[1], y1 = CIRCinside[2], n = 1)
    
     lines(dcirc)

     }
################

    ####   next get location of outside circle, circumcircle
    ###   the center is the intersection of perpendicular bi-sectors
    hdira = c(Ha[1]-P1[1], Ha[2]-P1[2])
    hdira = hdira/sqrt(sum(hdira^2))

    hdirb = c(Hb[1]-P2[1], Hb[2]-P2[2])
    hdirb = hdirb/sqrt(sum(hdirb^2))


    Q1 = c(U1[1]+ha*hdira[1], U1[2]+ha*hdira[2])
####   arrows(U1[1], U1[2], Q1[1], Q1[2])

    Q2 = c(U2[1]+hb*hdirb[1], U2[2]+hb*hdirb[2])
####  arrows(U2[1], U2[2], Q2[1], Q2[2])

    V2 = list(x=c(U2[1], Q2[1]), y=c(U2[2], Q2[2])) 
    V1 =  list(x=c(U1[1], Q1[1]), y=c(U1[2], Q1[2]) )
    CIRCUM = unlist(Sect2vex(V1, V2))
    
    if(add)  points(CIRCUM[1], CIRCUM[2])
  ##  circR = sqrt( (CIRCUM[1]-P3[1])^2 + (CIRCUM[2]-P3[2])^2  )
     if(add)
     {
       dcirc =  GEOmap::darc(rad = R, ang1 = 0, ang2 = 360, x1 =CIRCUM[1] , y1 =CIRCUM[2]  , n = 1)
    
     lines(dcirc)

     }
   
   ##### dcirc =  darc(rad = R, ang1 = 0, ang2 = 360, x1 = xp, y1 = yp, n = 1)

  #####  lines(dcirc)

  invisible(list(BI=CIRCinside, CIRCUM=CIRCUM,  IH=IH, CEN=CEN, r=r, R=R,AngBis=AngBis,
              H=c(ha, hb, hc), ALT=ALTs, M=c(ma, mb, mc), MED=MEDs, TEE=c(ta, tb, tc), Area=Area ))
  }
