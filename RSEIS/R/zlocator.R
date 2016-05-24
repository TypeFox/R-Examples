`zlocator` <-
function(COL=1, ID=FALSE, NUM=FALSE, YN=NULL, style=0)
  {
    ###   put a small vertical bar on the plot where you pick
    if(missing(COL)) { COL = 1 }
    if(missing(NUM)) { NUM = FALSE }
    if(missing(YN)) { YN = 1 }
    if(missing(style)) { style  = 0 }
    if(missing(ID)) { ID = FALSE }
    
    u = par("usr")
    # du = (u[4]-u[3])/(YN)
    du = 1/(YN)
    
    n = 0
    xsave = NULL
    ysave = NULL
    zloc = locator(1)
    while(length(zloc$x)>0)
      {
        xsave = c(xsave, zloc$x)
         ysave = c(ysave, zloc$y)

        ##  j = floor(YN*(zloc$y-u[3] )/(u[4]-u[3]))
        j = floor((zloc$y)/du)

        inout = (zloc$y>=u[3] & zloc$y<=u[4] &zloc$x>=u[1] &zloc$x<=u[2])
        
        y1 = j*du
        y2 = y1+du
        
#  print(paste(sep=' ', j, y1, y2))
          if(style==0)
            {
              abline(v=zloc$x[1], col=COL)
            }
          if(style==1)
            {
              segments(zloc$x, y1, zloc$x, y2, col=COL)
            }
          if(style==2)
            {
              abline(v=zloc$x[1], col=gray(0.88))
              segments(zloc$x, y1, zloc$x, y2, col=COL)
            }
        

       
        #  abline(v=zloc$x[1], col=COL)
        n = n+1
        if(NUM==TRUE)
          {
            text(zloc$x[1], u[4], labels=n, pos=3, xpd=TRUE)

          }

        ##  if(ID==TRUE & inout==TRUE)
        if(ID==TRUE)
          {
            
            alabs = format.default(zloc$x[1], digits=3)
            mtext(alabs,at=zloc$x[1], side=3, line=0, srt=45)

          }
        
        zloc = locator(1)
      }

    jj = YN-floor(ysave*YN)

    
    return(list(x=xsave, y=ysave, n=jj))
  }

