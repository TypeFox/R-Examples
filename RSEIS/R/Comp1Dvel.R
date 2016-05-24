`Comp1Dvel` <-
function(v, v2)
  {
    plot(c(v$vp,v$vs, v2$vp, v2$vs), c( -v$zp,-v$zs, -v2$zp,-v2$zs), type='n', xlab="Velocity, km/s", ylab="Depth, km")
     u = par('usr')
    lines(v$vp, -v$zp, type='s', col=4)
    lines(v$vs, -v$zs, type='s', col=3)

    lines(v2$vp, -v2$zp, type='s', col=4)
    lines(v2$vs, -v2$zs, type='s', col=3)

    grid()
    title(v$name)
  }

