`whichbutt` <-
function(v, buttons)
  {
  ##X## return which button was clicked   
    butts = rep(0, length(v$x))
    
    for(i in 1:length(butts))
      {
        
        KT = which(v$x[i]>buttons$x1 & v$x[i]<buttons$x2 & v$y[i]>buttons$y1 & v$y[i]<buttons$y2  )
        
        
        if(length(KT)>0)
          {
            
            butts[i] = KT
          }
      }
    
    return(butts)
    
  }

