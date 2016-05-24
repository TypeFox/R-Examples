`ymarginfo` <-
function(SIDE=1, s1=0.1, s2=0.8)
  {
    ### SIDE = 1, 3 (bottom, top)
    if(missing(SIDE)) { SIDE = 1 }
    if(missing(s1)) { s1=0.1 }
    if(missing(s2)) { s2=0.8 }
    
    u = par("usr")
    fin = par("fin")
    pin = par("pin")

    imarg = (fin[2]-pin[2])/2
    uinch =   (u[4]-u[3])/pin[2]

    if(SIDE==3)
      {
        if(par("ylog")==TRUE)
          {
            y1 = (10^(u[4]+ s1*imarg*uinch))
            y2 = (10^(u[4]+ s2*imarg*uinch))

          }
        else
          {
            y1 = (u[4]+ s1*imarg*uinch)
            y2 = (u[4]+ s2*imarg*uinch)
          }
      }
    else
      {
        if(par("ylog")==TRUE)
          {
            y1 = (10^(u[3]- s1*imarg*uinch))
            y2 = (10^(u[3]- s2*imarg*uinch))

          }
        else
          {
            y1 = ((u[3]- s1*imarg*uinch))
            y2 = ((u[3]- s2*imarg*uinch))
          }



      }

    
    
    return(c(y1, y2))
  }

