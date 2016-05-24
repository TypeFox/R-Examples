vecproj<-function(P1, P2)
  {
    ####   calculates the cos of the angle beteen two vectors
    if(is.list(P1))
      {
        P1 = unlist(P1)

      }
    if(is.list(P2))
      {
        P2 = unlist(P2)

      }

  A1 =   sqrt( sum(P1*P1) )
  A2 =   sqrt( sum(P2*P2) )

    if(A1<=0 | A2<=0 ) return(NULL)
    
    cang = sum(P1*P2)/(A1*A2)

    angrad = acos(cang)

    angdeg= angrad*180/pi

    d1  = cang*A2
    d2 =  cang*A1

    ###   angle in radians
    ###  dis1 = projected distance of P2  along P1
    ###  dis2 = projected distance of P1  along P2
    
    return(list(cang=cang, angrad=angrad, angdeg=angdeg, dis1=d1, dis2=d2 ))
  }
