sigconvGR <-
function(wigmat, wavepulse, dt)
  {
   
    d = dim(wigmat)
    Ntrace = d[2]

    NLEN=d[1]
    
    for(i in 1:Ntrace)
      {

        
        grlen = floor(.6/dt)
        fgr = 10
        tape = applytaper( rep(1, grlen), p = 0.2)
        tgr = seq(from=0, by=dt, length=grlen)
        siggr = tape*sin(2*pi*fgr*tgr)

        
        KL = floor( length(wavepulse)/2 )
        cx2 = convolve(wigmat[,i], wavepulse, conj = TRUE, type = c("open"))
        cx1 = cx2[(KL):((KL)+NLEN-1)]
        wigmat[,i] = cx1 

      }

    invisible(wigmat)

  }

