`vline` <-
function(x, per=1, COL=1, NUM=FALSE, LAB=1:length(x), lwd=0, lty=1)
  {
    if(missing(COL)) { col = gray(0.8) }
    if(missing(per)) { per=1  }
    if(missing(NUM)) { NUM = FALSE }
    if(missing(LAB)) { LAB = NULL }
    if(missing(lwd)) { lwd=1  }
    if(missing(lty)) { lty=1  }
    ##  if(missing()) { =1  }



    
    u = par("usr")
    n = length(x)
    dy = u[4]-u[3]
    if(per>=0)
      {
        y1 = u[3]
        y2 = y1+per*dy
      }
    else
      {
        y1 = u[4]
        y2 = y1+per*dy  
      }
    segments(x, rep(y1,n),  x, rep(y2, n), col=COL, lwd=lwd, lty=lty)
    
    if(NUM==TRUE | !is.null(LAB)  )
      {

        if(is.null(LAB) ) { LAB=1:length(x) }
        
        if(per>=0)
          {
             text(x, rep(y1,n), labels=LAB, pos=1, xpd=TRUE)
          }
        else
          {
             text(x, rep(y1,n), labels=LAB, pos=3, xpd=TRUE)
          }
      }
  }

