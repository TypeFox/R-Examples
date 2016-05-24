labelLine<-function(P1, P2, above=TRUE, dinch=.2, lab="text",
                    acode=3, alength=.06, aty=1, acol='black', bty=1, bcol='black', tcol='black', font=1, cex=1   )
  {

    if(missing(dinch)) { dinch=.2  }
    if(missing(lab)) { lab=""  }
    if(missing(acode)) { acode=3  }
    if(missing(alength)) { alength=.06  }
    if(missing(aty)) {  aty=1 }
    if(missing(acol)) {  acol='black'   }
    if(missing(above)) {  above=TRUE   }

    if(missing(bty)) {  bty=1 }
    if(missing(bcol)) {  bcol='black'   }

    if(missing(font)) {  font=1   }
    if(missing(cex)) {   cex=1   }
    if(missing(tcol)) {   tcol = 'black'   }

    
    k = 3
    
    v1 = c(P2[1]-P1[1] , P2[2]-P1[2])
    lenv1 = sqrt(sum(v1*v1))
    v1 = v1/lenv1

    

    perp1 = c(-v1[2] ,  v1[1]        )

    perp1 = perp1/sqrt(sum(  perp1*perp1  ))

    rotangle = atan2(v1[2], v1[1] )
    

       if(perp1[2]<0)
         {
           perp1 = -perp1

        rotangle = atan2(-v1[2], -v1[1] )
      }
    
    srt=rotangle*180/pi
 


    Wlab = strwidth(lab, units = "user", font=font, cex=cex)
    Hlab = strheight(lab, units = "user", font=font, cex=cex)
   
    
  ##    SENS1 = sign( AXB.prod(P2, P1)[3] )

  ## title(SENS1)

  ##  plot(c(P1[1], P2[1]), c(P1[2], P2[2] ), asp=1, type='n' )
   ##      segments(P1[1], P1[2], P2[1], P2[2])
    
    u = par("usr")
    pinch = par("pin")

    dx = dinch*(u[2]-u[1])/pinch[1]
    dy =  dinch*(u[4]-u[3])/pinch[2]

    textshift = 0.05 
    tdx =  textshift*(u[2]-u[1])/pinch[1]
    tdy =  textshift*(u[4]-u[3])/pinch[2]

    if(above==FALSE)
      {
        perp1 = -perp1
             
        tdx = Hlab+ textshift*(u[2]-u[1])/pinch[1]
        tdy = Hlab+ textshift*(u[4]-u[3])/pinch[2]
         
      }
    
    


    SH1 = c(P1[1]+dx*perp1[1], P1[2]+dy*perp1[2])
    SH2 =c(    P2[1]+dx*perp1[1] , P2[2]+dy*perp1[2]        ) 

    A1 = c( mean(c(SH1[1],P1[1])),  mean(c(SH1[2],P1[2])))
    A2 = c( mean(c(SH2[1],P2[1])),  mean(c(SH2[2],P2[2])))
    
   segments( P1[1], P1[2],  SH1[1] , SH1[2], lty=bty, col=bcol)
   segments( P2[1], P2[2],  SH2[1] , SH2[2], lty=bty, col=bcol)

    arrows(A1[1], A1[2], A2[1], A2[2], code=acode, length=alength, lty=aty, col=acol  )



   ##  text(mean( c(A1[1],A2[1])), mean( c(A1[2],A2[2]))  , labels=lab, pos=k)

    TX = mean( c(A1[1],A2[1]))+tdx*perp1[1]
    TY = mean( c(A1[2],A2[2]))+tdy*perp1[2]

    text(TX, TY , labels=lab, srt=srt, adj=c(.5, 0), font=font, cex=cex  )

  
    
    
  }
