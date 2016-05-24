`XSEC.drive` <-
  function(MOD, x1, y1, x2, y2 , zmax=100, COL=tomo.colors(100), LIM=NULL, STA=NULL, TOP=NULL, STDLAB=c("DONE", "QUIT"))
{

  if(missing(COL)) { COL=tomo.colors(100) }
  if(missing(LIM)) { LIM=NULL }
  if(missing(STA)) { STA=NULL }
  if(missing(TOP)) { TOP=NULL }
  if(missing(zmax)) { zmax = max(MOD$D) }

  if(missing(STDLAB))
    {
      STDLAB = c("DONE", "QUIT", "FRESH", "CLICKS", "AREA", "ASTATS", "LINE", "TOPO")
    }
  
  stdlab =  STDLAB
  
  BLABS = c(stdlab)
  NLABS = length(BLABS)
  NOLAB = NLABS +1000
  
  colabs = rep(1,length(BLABS))
  pchlabs = rep(4,length(BLABS))
  
  XZ = TOMOXSEC(MOD, x1, y1, x2, y2 , zmax=zmax, COL=COL, LIM=LIM, STA=STA, PLOT=FALSE)
  PLOT.TOMOXSEC(XZ, COL=COL)
  buttons = RPMG::rowBUTTONS(BLABS, col=colabs, pch=pchlabs)
  
  MAINdev = dev.cur()
  dev.set( MAINdev)
  u = par("usr")
  upar = par("usr")
  sloc = list(x=c(u[1],u[2]))

  iloc = locator(1, col=rgb(1,0.8, 0.8), type='p')
  zloc = iloc
  
  
  Nclick = length(zloc$x)
  
  if(is.null(iloc$x)) { return(0) }
  K = RPMG::whichbutt(iloc ,buttons)
  
  sloc = zloc

  
  PLOC = NULL
  
  while(TRUE)
    {
      if(K[Nclick] == match("DONE", BLABS, nomatch = NOLAB))
        {
           buttons = RPMG::rowBUTTONS(BLABS, col=rep(grey(.8), length(BLABS)), pch=rep("NULL", length(BLABS)))
          title("Return to Calling Program")
          
          break;
        }
      
      if(K[Nclick] == match("QUIT", BLABS, nomatch = NOLAB))
        {

           buttons = RPMG::rowBUTTONS(BLABS, col=rep(grey(.8), length(BLABS)), pch=rep("NULL", length(BLABS)))
          title("Return to Calling Program")
          
          return(NULL)
        }


      if(iloc$x<upar[1] & (iloc$y>upar[3] & iloc$y<upar[4]) )
        {
          

          PLOT.TOMOXSEC(XZ, COL=COL)
          upar = par("usr")
          buttons = RPMG::rowBUTTONS(BLABS, col=colabs, pch=pchlabs)
          zloc = list(x=NULL, y=NULL)
          
        }
      
      
      if(K[Nclick] == match("FRESH", BLABS, nomatch = NOLAB))
        {
          
          PLOT.TOMOXSEC(XZ, COL=COL)
          upar = par("usr")
          buttons = RPMG::rowBUTTONS(BLABS, col=colabs, pch=pchlabs)
          zloc = list(x=NULL, y=NULL)
          
        }
      if(K[Nclick] == match("TOPO", BLABS, nomatch = NOLAB))
        {
          lines(TOP$x, TOP$z)
          zloc = list(x=NULL, y=NULL)

        }

      iloc = locator(1, col=rgb(1,0.8, 0.8), type='p')
      zloc  = list(x=c(zloc$x,iloc$x), y=c(zloc$y, iloc$y))
      Nclick = length(iloc$x)
      
      if(is.null(iloc$x)) { return(zloc) }
      K =  RPMG::whichbutt(iloc ,buttons)
      
      
      
    }

}

