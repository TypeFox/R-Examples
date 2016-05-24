`winmark` <-
  function(a1, a2, side=1, bar=NULL, leg=NULL, col=col, lwd=1, lty=1, arrows=FALSE, alen=.1, leglen=0.15, LEGON=3, BARON=TRUE )
{
  if(missing(col)) { col=4 }
   if(missing(lwd)) { lwd=1 }
   if(missing(lty)) { lty=1 }
   if(missing(arrows)) { arrows=FALSE }
   if(missing(alen)) { alen=0.125  }
   if(missing(leglen )) { leglen=0.125  }
  if(missing(BARON)) { BARON = TRUE }
  if(missing(LEGON)){ LEGON= 3 }
    
  if(missing(side)) { side=1 }

  if(side<1 | side>4) side = 1

  margpct   = 0.05
  legpct =  leglen
  
  if(missing(leg) | missing(bar)) {
    u = par("usr")
    if(side==3)
      {
        leg = u[3]+(1-legpct) *(u[4]-u[3])
        bar = u[3]+ (1-margpct) *(u[4]-u[3])
        
      }
     if(side==1)
      {
        
        bar = u[3]+margpct*(u[4]-u[3])
        leg  = u[3]+legpct*(u[4]-u[3])
        
      }

    if(side==4)
      {
        leg = u[1]+(1-legpct)*(u[2]-u[1])
        bar = u[1]+(1-margpct)*(u[2]-u[1])
        
      }
     if(side==2)
      {
        
        bar = u[1]+(margpct)*(u[2]-u[1])
        leg = u[1]+(legpct)*(u[2]-u[1])
        
      }
    
  }

###  bar[ is.na(bar) ] = 1
###  leg[ is.na(leg) ] = 1
  

  if(is.null(a1)==TRUE )
    {
      print("missing a1 in winmark")
      return(0)
    }

  if(is.null(a2)==TRUE )
    {
      print("missing a2 in winmark")
      return(0)
    }

 ###### print(bar)
 ###### print(leg)

###############  take careof conflicts in bar, leg and side
######   if bar>leg but side = 1 or 2 this is a conflict
######   if bar<leg  but side = 3,4   this is a conflict
    
####  provision of bar and leg over-rules side

  if(side==1 &  all(bar>leg) )
    {
      side =3
    }
  if(side==2 &  all(bar>leg) )
    {
      side =4
    }
  if(side==3 &  all(bar<leg) )
    {
      side =1
    }
  if(side==4 &  all(bar<leg) )
    {
      side =2
    }
  
  #########  if side == 1 then bar is on bar and legs point down

  if(side ==  1 | side ==3)
    {

      if(LEGON==1 | LEGON==3) segments(a1, leg, a1, bar, col=col, lwd=lwd, lty=lty, xpd=TRUE)
      if(BARON)               segments(a1, bar, a2, bar, col=col, lwd=lwd, lty=lty, xpd=TRUE)
      if(LEGON==2 | LEGON==3) segments(a2, bar, a2, leg, col=col, lwd=lwd, lty=lty, xpd=TRUE)

      if(arrows)
        {
       if(LEGON==1 | LEGON==3)   arrows(a1, bar, a1, leg, col=col, lwd=lwd, length=alen)
       if(LEGON==2 | LEGON==3)    arrows(a2, bar, a2, leg, col=col, lwd=lwd, length=alen)
          
        }
    }
  else
    {


      if(LEGON==1 | LEGON==3) segments( leg, a1, bar, a1, col=col, lwd=lwd, lty=lty, xpd=TRUE)
      if(BARON)  segments( bar, a2, bar, a1, col=col, lwd=lwd, lty=lty, xpd=TRUE)
      if(LEGON==2 | LEGON==3)segments( bar, a2, leg, a2, col=col, lwd=lwd, lty=lty, xpd=TRUE)

      if(arrows)
        {
        if(LEGON==1 | LEGON==3)   arrows( bar, a1, leg, a1, col=col, lwd=lwd, length=alen)
       if(LEGON==2 | LEGON==3)    arrows( bar, a2, leg, a2, col=col, lwd=lwd, length=alen)
          
        }

    }

}

