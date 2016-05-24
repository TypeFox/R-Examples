xsecmanyfoc<-function(MEK, focsiz=0.04, theta=NULL, foccol=NULL, UP=TRUE, LEG = FALSE, DOBAR=FALSE)
  {
    if(missing(theta)) { theta = NULL }
    
    if(missing(UP)) UP = TRUE
  
    
    ############################  focal mech have this list structure

#########    MEKS = list(lon=0, lat=0, str1=0, dip1=0, rake1=0, lat=0, lon=0, dep=0, name="", Elat=0, Elon=0)
    ###  Elat and Elon are the exploded lat-lon (a line will be drawn from the
    ########    epicenter to the plotting position
    
   
    if(missing(foccol)) foccol=NULL
    if(missing(LEG))   LEG = FALSE
    if(missing(DOBAR))  DOBAR=FALSE

    
    if(missing(focsiz)) { focsiz = 0.04 }
  
if(length(MEK$Elat)<1 | (length(MEK$Elon)<1))
  {
    MEK$Elat = rep(NA, length(MEK$lat))
    MEK$Elon = rep(NA, length(MEK$lat))

  }

   u = par("usr")
  
    
    userp = par("usr")
    useri = par("pin")

    ratx = (userp[2]-userp[1])/useri[1]
    raty=  (userp[4]-userp[3])/useri[2]

  
   
    focsizx = focsiz*ratx
    focsizy = focsiz*raty
    
    
    XYplode = list(x=MEK$Elat, y=MEK$Elon)

    tem1 = list(x=MEK$x, y=MEK$y)
    
    tem2 = XYplode


    ww = which(!is.na(MEK$x))


    
    if(  length(ww) < 1 ) {   return()  }

      ##############################  loop through focal mechs and plot each one
    for(IW in 1:length(ww))
      {

        i = ww[IW]
        

        lind= foc.icolor(MEK$rake1[i])
        focal.col = foc.color(lind, pal=1)

        
        a1  =  SDRfoc(MEK$str1[i], MEK$dip1[i],MEK$rake1[i] , u=UP, ALIM=c(-1,-1, +1, +1), PLOT=FALSE)

        if(!is.null(theta))
          {
            
            b1 = Rotfocphi(theta, a1$M$uaz, a1$M$ud, a1$M$vaz,
              a1$M$vd,  a1$M$az1,  a1$M$d1,  a1$M$az2,  a1$M$d2,
              a1$M$paz,  a1$M$pd,  a1$M$taz,  a1$M$td)


 
            MEC=list(UP=a1$UP, LIM=a1$LIM, P=b1$P,  az1=b1$A2$az-90  , 
              dip1=b1$A2$dip , az2=b1$A1$az-90, dip2=b1$A1$dip )

            MEC$sense = 0
            
            pax = focpoint(MEC$P$az, MEC$P$dip,  lab="P", UP=MEC$UP, PLOT=FALSE)
  
            PLS = polyfoc(MEC$az1, MEC$dip1, MEC$az2, MEC$dip2, UP=MEC$UP, PLOT = FALSE)
            kin = inout(cbind(pax$x, pax$y) ,cbind(PLS$Px, y =PLS$Py), bound=TRUE)

            if(kin==0) MEC$sense = 1
            
          }
        else
          {
            MEC = a1

          }
        

       if( is.na(tem2$x[i]) |  is.na(tem2$y[i]))
          {
            
            justfocXY(MEC, x=tem1$x[i],  y= tem1$y[i],  size = c(focsizx, focsizy), fcol =focal.col , fcolback = "white", xpd = TRUE)
          }
          else
          {
           ##  print(paste("Exploding, i=", i))
            segments(tem1$x[i], tem1$y[i], tem2$x[i], tem2$y[i])
            justfocXY(MEC, x=tem2$x[i],  y= tem2$y[i],  size = c(focsizx, focsizy), fcol =focal.col , fcolback = "white", xpd = TRUE)


            if(DOBAR)
              {

                StrikeDip( tem1$x[i], tem1$y[i] , MEC  ,focsizx, addDIP=FALSE, col="black", lwd=3 )
              }


          }
      }



    if(LEG)
      {
        fleg = 1:7
        flegc = foc.color(fleg, pal=1)
        flab = focleg(fleg)

        legend("topleft", legend=flab, col=flegc, pch=19, bg="white" )
      }


    
  }

