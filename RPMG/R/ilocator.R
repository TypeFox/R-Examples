`ilocator` <-
function(N=1, COL=1, NUM=FALSE, YN=NULL, style=0)
  {
    ###   put a small vertical bar on the plot where you pick
    if(missing(COL)) { COL = 1 }
    if(missing(NUM)) { NUM = FALSE }
    if(missing(YN)) { YN = 1 }
    if(missing(style)) { style  = 0 }
    if(missing(N)) { N=1 }

    u = par("usr")
    # du = (u[4]-u[3])/(YN)
    du = 1/(YN)
    
    n = 0
    xsave = NULL
    ysave = NULL
    ## zloc = locator(1)
    for(iK in 1:(N))
      {
        zloc = locator(1)
        if(length(zloc$x)<1) return(NULL) 
        xsave = c(xsave, zloc$x)
         ysave = c(ysave, zloc$y)

        ##  j = floor(YN*(zloc$y-u[3] )/(u[4]-u[3]))
        j = floor((zloc$y)/du)
       
        y1 = j*du
        y2 = y1+du
        
#  print(paste(sep=' ', j, y1, y2))

        if(style==(-1))
          {
            points(zloc$x, zloc$y,  pch=23, col=COL)
          }     
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
          if(style==3)
            {
              abline(v=zloc$x[1], col=gray(0.88))
              segments(zloc$x, y1, zloc$x, y2, col=COL)
              segments(u[1], zloc$y, u[2], zloc$y, col=COL)
              
            }
        

       
        #  abline(v=zloc$x[1], col=COL)
        n = n+1
        if(NUM==TRUE)
          {
            if(style==(-1))
              {
                text(zloc$x[1], zloc$y[1], labels=n, pos=3, xpd=TRUE)
              }
            else
              {
                text(zloc$x[1], u[4], labels=n, pos=3, xpd=TRUE)
          }

          }
        
      }

    jj = YN-floor(ysave*YN)

    
    return(list(x=xsave, y=ysave, n=jj))
  }

