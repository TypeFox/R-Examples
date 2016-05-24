randFRY<-function(N=200, LIM=c(0,0, 100, 100) , rlen=5, PLOT=FALSE)
  {
    if(missing(N)) { N=200 }
    if(missing(LIM))  LIM=c(0,0, 100, 100)
    if(missing(rlen))  rlen=5
    if(missing(PLOT)) PLOT=FALSE

    
    x = vector()
    y = vector()

    px = runif(1, LIM[1], LIM[3])
    py = runif(1, LIM[2], LIM[4])

    x[1] = px
    y[1] = py

    k = 1
    
    while(k<N)
      {
        px = runif(1, LIM[1], LIM[3])
        py = runif(1, LIM[2], LIM[4] )

        r = sqrt( (px-x)^2 + (py-y)^2)
        if(all(r>rlen))
          {
            k = k+1
            x[k] = px
            y[k] = py


          }

      }


    if(PLOT) plot(x,y, asp=1, pch=".", cex=3)

    return(list(x=x, y=y) )
  }
