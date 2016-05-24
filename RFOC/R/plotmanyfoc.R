plotmanyfoc<-function(MEK, PROJ, focsiz=0.04, foccol=NULL, UP=TRUE , PMAT=NULL, LEG = FALSE, DOBAR=FALSE)
  {

    ############################  focal mech have this list structure

#########    MEKS = list(lon=0, lat=0, str1=0, dip1=0, rake1=0, lat=0, lon=0, dep=0, name="", Elat=0, Elon=0)
    ###  Elat and Elon are the exploded lat-lon (a line will be drawn from the
    ########    epicenter to the plotting position
    
    if(missing(PMAT)) PMAT=NULL
    if(missing(foccol)) foccol=NULL
    if(missing(UP)) UP=TRUE
    
    if(missing(LEG))   LEG = FALSE
    if(missing(DOBAR))  DOBAR=FALSE

    
    if(missing(focsiz)) { focsiz = 0.04 }
  
if(length(MEK$Elat)<1 | (length(MEK$Elon)<1))
  {
    MEK$Elat = rep(NA, length(MEK$lat))
    MEK$Elon = rep(NA, length(MEK$lat))

  }

    ## require(GEOmap)
    
    u = par("usr")
    
    focsiz = focsiz* (u[2]-u[1])
    
    
    EQxy = GEOmap::GLOB.XY(MEK$lat, MEK$lon, PROJ)

    XYplode = GEOmap::GLOB.XY(MEK$Elat, MEK$Elon, PROJ)

    tem1 = EQxy
    tem2 = XYplode
    if(!is.null(PMAT))
      {
        tem1 = trans3d(EQxy$x, EQxy$y, rep(0, length(EQxy$y)), PMAT)
        tem2 = trans3d(XYplode$x, XYplode$y, rep(0, length(XYplode$y)), PMAT)
      }
    
    for(i in 1:length(MEK$str1))
      {

       lind= foc.icolor(MEK$rake1[i])
       focal.col = foc.color(lind, pal=1)

        
        MEC =  SDRfoc(MEK$str1[i], MEK$dip1[i],MEK$rake1[i] , u=UP, ALIM=c(-1,-1, +1, +1), PLOT=FALSE)

       if( is.na(tem2$x[i]) |  is.na(tem2$y[i]))
          {
            
            justfocXY(MEC, x=tem1$x[i],  y= tem1$y[i],  size = c(focsiz, focsiz), fcol =focal.col , fcolback = "white", xpd = TRUE)
          }
          else
          {
           ##  print(paste("Exploding, i=", i))
            segments(tem1$x[i], tem1$y[i], tem2$x[i], tem2$y[i])
            justfocXY(MEC, x=tem2$x[i],  y= tem2$y[i],  size = c(focsiz, focsiz), fcol =focal.col , fcolback = "white", xpd = TRUE)


            if(DOBAR)
              {

                StrikeDip( tem1$x[i], tem1$y[i] , MEC  ,focsiz, addDIP=FALSE, col="black", lwd=3 )
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

