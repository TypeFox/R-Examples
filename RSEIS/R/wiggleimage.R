wiggleimage <-
function(Arot, dt=1, dx=1, col='black' )
  {
    if(missing(dt)) { dt = 1 }
    if(missing(dx)) { dx = 1 }
    if(missing(col)) { col='black' }
    

    d = dim(Arot)


    if(length(dx)>1)
      {
        Ax = dx
        dx = dx[2]-dx[1]
      }
    else
      {
        Ax = seq(from=0, by=dx, length=(d[2]))

      }

    if(length(dt)>1)
      {
        Atim = dt
        dt = Atim[2]-Atim[1]
        
      }
    else
      {
    
    Atim = seq(from =0, by=dt,  length=(d[1]))
  }

    plot(range(c(0,Ax)), range(c(0, Atim)), type='n', xlab="km", ylab="time")
   ##  abline(h = seq(from=0, to=min(Atim) , by=(-0.1)*dt), col=grey(0.8))
    grid()


    segments(Ax, rep(0, times=length(Ax)), Ax, rep(max(Atim) , times=length(Ax)), col=grey(.8) )

    abline(h = 0, v=0)
    axis(3, at=Ax, pos=max(Atim))


    for(i in 1:length(Ax))
      {
        lines(Ax[i]+Arot[ ,i ], Atim, col=col)
      }

  }

