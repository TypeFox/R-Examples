NinePointCircle<-function (P1, P2 = c(0, 1), P3 = c(1, 0), add = FALSE, SHOW=TRUE) 
{

  if(missing(add)){ add=FALSE }
  if(missing(SHOW)) { SHOW=FALSE } 
  if (missing(P2) | missing(P3)) {
    if (is.matrix(P1)) {
      P2 = P1[2, ]
      P3 = P1[3, ]
      P1 = P1[1, ]
    }
    if (is.list(P1)) {
      P2 = c(P1$x[2], P1$y[2])
      P3 = c(P1$x[3], P1$y[3])
      P1 = c(P1$x[1], P1$y[1])
    }
  }



  arrows1=FALSE

 
  if(FALSE)
    {

      source("/Users/lees/RPAX/geophys/NinePointCircle.R")

      P1 = 10*runif(2)
      P2 =  10*runif(2)
      P3 =  10*runif(2)

      TRI =  NinePointCircle(P1, P2, P3, add=TRUE, SHOW=TRUE)

    }


  theX = c(P1[1],P2[1], P3[1])
  theY = c(P1[2],P2[2], P3[2])

  TRI =  TriangleInfo(P1, P2, P3, add=FALSE)

     #####   get the medians and find their circumscribed circle=nine point circle

  NineCent =  TriangleInfo(TRI$MED, add=FALSE)

  #######   altitude vectors
  H1 =  c((TRI$ALT[1,1]-P1[1]), (TRI$ALT[1,2]-P1[2])  )
  H2 =  c((TRI$ALT[2,1]-P2[1]), (TRI$ALT[2,2]-P2[2])  )
  H3 =  c((TRI$ALT[3,1]-P3[1]), (TRI$ALT[3,2]-P3[2])  )

  ###  S-points
  I2 = Sect2vex(list(x=c(TRI$ALT[1,1],P1[1])   , y=c(TRI$ALT[1,2],P1[2]) ),
    list(x=c(TRI$ALT[2,1],P2[1]), y=c(TRI$ALT[2,2],P2[2])))

     J =  c( I2$x+ 0.5*(P1[1]-I2$x)  , I2$y+ 0.5*(P1[2]-I2$y))
     K =  c( I2$x+ 0.5*(P2[1]-I2$x)  , I2$y+ 0.5*(P2[2]-I2$y))
     L =  c( I2$x+ 0.5*(P3[1]-I2$x)  , I2$y+ 0.5*(P3[2]-I2$y))


  
 if(SHOW==TRUE)
   {
     plot(c(theX,TRI$ALT[,1],I2$x    ), c(theY, TRI$ALT[,2],I2$y ) , type='n', asp=1, ann=FALSE)
   }


  if(add==TRUE)
    {
     
     points(theX , theY)
     text(theX , theY, labels=c("A", "B", "C") , pos=3)
     lines(c(theX, theX[1]), c(theY, theY[1]) )


     if(arrows1)
       {
         arrows(P1[1],P1[2], TRI$ALT[1,1],  TRI$ALT[1,2], length=.1)
         arrows(P2[1],P2[2], TRI$ALT[2,1],  TRI$ALT[2,2], length=.1)
         arrows(P3[1],P3[2], TRI$ALT[3,1],  TRI$ALT[3,2], length=.1)
       }
     
     points( TRI$MED[,1]  ,  TRI$MED[,2]  )

     text(TRI$MED[,1]  ,  TRI$MED[,2], labels=c("D", "E", "F"), pos=3)
     points(NineCent$CIRCUM[1], NineCent$CIRCUM[2], col='red')

     points(TRI$ALT[,1]  ,  TRI$ALT[,2], col='blue')

     text(TRI$ALT[,1]  ,  TRI$ALT[,2], labels=c("G", "H", "I"), pos=3)

     dcirc = GEOmap::darc(rad = NineCent$R, ang1 = 0, ang2 = 360, x1 = NineCent$CIRCUM[1], 
       y1 = NineCent$CIRCUM[2], n = 1)
     lines(dcirc)

     points(I2$x, I2$y, pch=3)
     text(I2$x, I2$y,  labels="s", pos=1)

     segments(I2$x, I2$y, P1[1],P1[2], col="gold")
     segments(I2$x, I2$y, P2[1],P2[2], col="gold")
     segments(I2$x, I2$y, P3[1],P3[2], col="gold")


     ##  j1= sqrt((P1[1]-I2$x)^2+(P1[2]-I2$y)^2)

     points(J[1], J[2])
     points(K[1], K[2])
     points(L[1], L[2])
     text(J[1], J[2], labels="J", pos=4)
     text(K[1], K[2], labels="K", pos=4)
     text(L[1], L[2], labels="L", pos=4)

     
     points(NineCent$CIRCUM[1], NineCent$CIRCUM[2], col='red', pch=3)
   }




 return(
         list(
              A=P1, B=P2, C=P3,
              D=c(TRI$MED[1,1]  ,  TRI$MED[1,2]),
              E=c(TRI$MED[2,1]  ,  TRI$MED[2,2]),
              F=c(TRI$MED[3,1]  ,  TRI$MED[3,2]),
              G=c(TRI$ALT[1,1]  ,  TRI$ALT[1,2]),
             H=c(TRI$ALT[2,1]  ,  TRI$ALT[2,2]),
             I=c(TRI$ALT[3,1]  ,  TRI$ALT[3,2]),
              J=J, K=K, L=L, S=I2, CEN=NineCent$CIRCUM,   R=NineCent$R))


         
              
}



