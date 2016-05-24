`VEC2MOD` <-
  function(VEC)
  {


    MOD = list()
    MOD$x = attr(VEC,'x')
    MOD$y = attr(VEC,'y')
    MOD$D = attr(VEC,'D')

    nx = length(MOD$x)
    ny = length(MOD$y)
    nz = length(MOD$D)
    lentop = nx*ny

    for(i in 1:nz)
      {
        j1=lentop*(i-1)+1
        j2=(j1-1)+lentop

        lilvec = VEC[j1:j2]
        lilvec[lilvec==0] = NA

        MOD$MOD[[i]] = matrix(lilvec, ncol=nx, nrow=ny)


      }
    invisible(MOD)

    
  }

