`Get1Dvel` <-
function(infile, PLOT=TRUE)
  {
    ###  read in a UW format velocity file
    if(missing(PLOT)) { PLOT=TRUE }

    if(!is.list(infile) & is.character(infile) )
      {   
        
        
        av = scan(file=infile, nlines = 2, what="" , sep="\n", quiet =TRUE)
        
        v =  scan(file=infile, skip=2, list(zp=0, vp=0, ep=0, zs=0, vs=0, es=0), quiet =TRUE)
        v$name = infile
        
        v$descriptor = av
      }
    else
      {
        ###########  velocity file  is already in memory.....
        v = infile
      }
    
    
                                        #
    if(PLOT)
      {
        plot(c(v$vp,v$vs), c( -v$zp,-v$zs), type='n', xlab="Velocity, km/s", ylab="Depth, km") 
        lines(v$vp, -v$zp, type='s', col=4)
        lines(v$vs, -v$zs, type='s', col=3)
        title(v$name)
        u = par('usr')
        LEG = legend("bottomleft"  , c("Vp", "Vs"), lwd=2, col=c(4,3), plot=TRUE  )
       ##### jlegend( u[1], u[3]+LEG$rect$h, c("Vp", "Vs"), lwd=2, col=c(4,3), plot=TRUE  )
      }
   
    return(v)
  }

