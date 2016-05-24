dowiggles <-
function(AMAT, dt, dx )
  {
    d = dim(AMAT)

    if(length(dx)>1)
      {
        x = dx
        dx = x[2]-x[1]
      }
    else
      {
        
        x = seq(from=0, by=dx, length=d[2])
      }


    
    x = seq(from=0, by=dx, length=d[2])

    PM =  matsquiggle(AMAT, dt, dist = x, thick = .09, FLIP = TRUE, filcol=NA , tracecol="black", add=FALSE, PLOT=FALSE)

    tmax = d[1]*dt
    axis(1, pos=0)
    axis(2)
    abline(h = 0, v=0)
    axis(3, at=x, pos=tmax)
    ##  abline(h = seq(from=0, to=PM$rx[2], by=0.1), col=grey(0.8))

    PM =  matsquiggle(AMAT, dt, dist = x, thick = .09, FLIP = TRUE, filcol=NA, tracecol="black", add=TRUE, PLOT=TRUE)


  }

