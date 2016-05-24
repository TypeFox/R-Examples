`Comp1Dvels` <-
function(INV, depth=1:50 )
  {
    #@####  compare several 1D velocity models
    if(missing(depth)) { depth=NULL }
    
    raj = vector()
    for(i in 1:length(INV))
      {
         v = get(INV[i])
         r1 = range(v$zp)
         r2 = range(v$vp)
         r3 = range(v$zs)
         r4 = range(v$vs)
         
          raj = rbind(raj, c(r1,r2,r3,r4 ) )
         
       }
    
    ###      plot(c(v$vp,v$vs, v2$vp, v2$vs), c( -v$zp,-v$zs, -v2$zp,-v2$zs), type='n', xlab="Velocity, km/s", ylab="Depth, km")

    D = -range(raj[,c(1,2,5,6) ])
    if(!is.null(depth)) {  D = -depth }
    
    plot(  range(raj[,c(3,4,7,8) ]), D, type='n', xlab="Velocity, km/s", ylab="Depth, km")
    
    
    u = par('usr')
    for(i in 1:length(INV))
      {
        v = get(INV[i])

      zp = c(v$zp, v$zp[length(v$zp)]+10)
      zs = c(v$zs, v$zs[length(v$zs)]+10)
        vp = c(v$vp[1], v$vp)
        vs = c(v$vs[1], v$vs)
        lines(vp, -zp, type='s', col=i+1, lwd=2)
        lines(vs, -zs, type='s', lty=2, col=i+1, lwd=2)
        
      }
    
    
    grid()
    legend("topright", legend =INV, lty=1, lwd=2, col=2:(length(INV)+1))

    
  }

