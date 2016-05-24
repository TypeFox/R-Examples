`see.pal` <-
function(col)
  {
  
    u = par("usr")
    f = par("pin")

    raty = (u[4]-u[3])/f[2]

    
    dy = (u[4]-u[3])*.05


    
    dx = (u[2]-u[1])*.3

    LU=list(x=c(u[1]+dx*0.1, u[1]+dx*0.1+dx), y = c(u[3]-dy*1.5, u[3]-dy*1.5-dy))
  
    
    i <- seq(along = col)
    BX = (LU$x[2]-LU$x[1])/length(i)
    x1 =LU$x[1]+(i-1)*BX
    x2 = x1+BX
    y1 = LU$y[1]
    y2 =  LU$y[2]
      
    rect(x1,y1,x2,y2,  col=col, xpd = TRUE, border=NA)
    
    invisible(c(x1,y1,x2,y2))
  }

