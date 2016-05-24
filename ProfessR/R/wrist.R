`wrist` <-
function(DB)
  {
    abline(v=DB$divs, col='blue')
    
    ex = (DB$divs[2:(length(DB$divs))] + DB$divs[1:(length(DB$divs)-1)])/2
    u =  par("usr")
    text( ex,  rep(u[4], length(ex)) , labels=c("E", "D", "C", "B", "A"), pos=1 , col=rgb(1, .7, .7)   )
    

  }

