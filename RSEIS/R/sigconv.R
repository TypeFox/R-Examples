sigconv <-
function(wigmat, wavepulse)
  {
   
    d = dim(wigmat)
    Ntrace = d[2]

    NLEN=d[1]
    KL = floor( length(wavepulse)/2 )
    for(i in 1:Ntrace)
      {

        cx2 = convolve(wigmat[,i], wavepulse, conj = TRUE, type = c("open"))
        cx1 = cx2[(KL):((KL)+NLEN-1)]
        wigmat[,i] = cx1 

      }

    invisible(wigmat)

  }

