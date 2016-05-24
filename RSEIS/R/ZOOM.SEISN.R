`ZOOM.SEISN` <-
function(GH, sel=1:length(GH$dt), WIN=NULL)
{
  if(missing(WIN)) { WIN = NULL }
  if(missing(sel)) { sel = 1:length(GH$dt)}

  if(is.logical(sel)) { sel= which(sel)  }
  
  labs = c("STOP", "zoom out", "zoom in", "restore", "saveWIN")
  
  YN = PLOT.SEISN(GH, WIN=WIN, dt=GH$dt, sel=sel , notes=GH$KNOTES[sel])
  u = par("usr")
  sloc = list(x=c(u[1],u[2]))
  #### ftime = Zdate(GH$info, sel[1],0)
  #### mtext( ftime, side = 3, at = 0, line=0.5, adj=0)
  
 ####  title("LEFT 0 Click = done; 1 Click=replot;   2 Click=zoom")
  buttons = RPMG::rowBUTTONS(labs)
  
  ####  NV = LabelBAR(labs)
  zloc = plocator(COL=rgb(1,0.8, 0.8))
  Nclick = length(zloc$x)
  if(is.null(zloc$x)) { return(NULL) }
  K = RPMG::whichbutt(zloc ,buttons)

  sloc = zloc

  while(Nclick>0)
    {

      if(K[Nclick] == 1)
        {
          break;
        }

      
      if(Nclick==1 & K[Nclick]==0)
        {
          WIN = NULL
          YN = PLOT.SEISN(GH, WIN=WIN, dt=GH$dt, sel=sel , notes=GH$KNOTES[sel])
            
          u = par("usr")
          
          sloc = list(x=c(u[1],u[2]))
          
        ####  ftime = Zdate(GH$info, sel[1], 0)
        ####  mtext( ftime, side = 3, at = 0, line=0.5, adj=0)
        ####  next;
          
        }

      if(K[Nclick]==4)
        {
          WIN = NULL
          YN = PLOT.SEISN(GH, WIN=WIN, dt=GH$dt, sel=sel , notes=GH$KNOTES[sel])
          u = par("usr")
          L = length(sloc$x)
          if(L>1)
            {
              abline(v=sloc$x[c(L-1,L)], col=gray(0.8), lty=2)
            }
          sloc = list(x=c(u[1],u[2]))
          
        }

      if(Nclick>1 & K[Nclick]==0)
        {

          WIN  = zloc$x[c( Nclick-1, Nclick)]
          YN = PLOT.SEISN(GH, WIN=WIN, dt=GH$dt, sel=sel , notes=GH$KNOTES[sel])
         #### ftime = Zdate(GH$info, sel[1], WIN[1])
         ####  mtext( ftime, side = 3, at = WIN[1], line=0.5, adj=0)
          sloc = zloc
        }
      
      if(K[Nclick]==2)
        {
          u = par("usr")
          DX = (u[2]-u[1])*0.3
          zloc = list(x= c(u[1]-DX, u[2]+DX))
          WIN  = zloc$x
          YN = PLOT.SEISN(GH, WIN=WIN, dt=GH$dt, sel=sel , notes=GH$KNOTES[sel])
         #### ftime = Zdate(GH$info, sel[1], WIN[1])
         ####  mtext( ftime, side = 3, at = WIN[1], line=0.5, adj=0)

          sloc = zloc
        }
       if(K[Nclick]==3)
        {
          u = par("usr")
          DX = (u[2]-u[1])*0.3
          zloc = list(x= c(u[1]+DX, u[2]-DX))
          WIN  = zloc$x
         
          YN = PLOT.SEISN(GH, WIN=WIN, dt=GH$dt, sel=sel , notes=GH$KNOTES[sel])
         #### ftime = Zdate(GH$info, sel[1], WIN[1])
         ####  mtext( ftime, side = 3, at = WIN[1], line=0.5, adj=0)

          sloc = zloc
        }
       if(K[Nclick]==5)
        {
          
          print(paste(sep=" " , "WIN=",sloc$x))
          

          
        }
          
      buttons = RPMG::rowBUTTONS(labs)
     ###  NV = LabelBAR(labs)
      zloc = plocator(COL=rgb(1,0.8, 0.8))
      Nclick = length(zloc$x)
      if(is.null(zloc$x)) { return(sloc) }
      K =  RPMG::whichbutt(zloc ,buttons)
      ### K = ValBAR(NV, zloc)
     ###  print(paste(sep=" ", "K=",K))
      
    }
  return(sloc)
}

