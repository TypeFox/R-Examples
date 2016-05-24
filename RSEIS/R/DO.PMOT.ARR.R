`DO.PMOT.ARR` <-
function(E, N)
  {
    len = floor( (length(N)-4) / 4  )
    for(jarr in 1:4)
      {
        karr = (jarr)*len
        harr = hypot(E[karr], N[karr], E[karr+1], N[karr+1])
       #### print(c(jarr, karr, harr))
        
       ######  print(paste(sep=' ', jarr, karr, harr))
        if(!is.na(harr) & harr>.001)
          {
            arrows(E[karr], N[karr], E[karr+1], N[karr+1], length=.06)
          }
      }

  }

