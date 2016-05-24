`piczoeppritz` <-
function(LL=list(x=c(0,1), y=c(0,1)) , chincw="P")
  {
    if(missing(LL)) { LL = locator(2) }
    if(missing(chincw)) { chincw="P" }

     zoepcols = c("red", "green" , "blue", "purple")

    mx = mean(LL$x)
    my =  mean( LL$y)
    rect(LL$x[1], LL$y[1], LL$x[2], LL$y[2],  col=NULL, border=grey(0.85) )
    
    segments(LL$x[1], my ,  LL$x[2], my, lty=2  )  
    segments(mx, LL$y[1],  mx, LL$y[2], lty=2   )
    
    arrows(LL$x[1]+0.25*(diff(LL$x))  , LL$y[2], mx, my, length = 0.1 )
    text(LL$x[1]+0.25*(diff(LL$x))  , LL$y[2], labels=chincw, pos=3)
    
    
    arrows(mx, my, LL$x[1]+0.65*(diff(LL$x)), LL$y[2], length = 0.1, col=zoepcols[2])
    arrows(mx, my, LL$x[1]+0.85*(diff(LL$x)), LL$y[2], length = 0.1, col=zoepcols[1])

    text( LL$x[1]+0.65*(diff(LL$x)), LL$y[2]   , labels="S", pos=3, col=zoepcols[2])
      text(   LL$x[1]+0.85*(diff(LL$x)), LL$y[2]  , labels="P", pos=3, col=zoepcols[1])
    

    arrows(mx, my, LL$x[1]+0.65*(diff(LL$x)), LL$y[1], length = 0.1, col=zoepcols[4])
    arrows(mx, my, LL$x[1]+0.85*(diff(LL$x)), LL$y[1], length = 0.1, col=zoepcols[3])

 text(LL$x[1]+0.65*(diff(LL$x)), LL$y[1],  labels="S", pos=1, col=zoepcols[4])
 text( LL$x[1]+0.85*(diff(LL$x)), LL$y[1],   labels="P", pos=1, col=zoepcols[3])

    
  }

