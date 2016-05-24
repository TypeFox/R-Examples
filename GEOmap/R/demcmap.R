`demcmap`<-function(ZTOPO, n=100, ccol=NULL)
  {
    if(missing(n)) n = 100

    ####   ZTOPO is a vector of elevations - can be length 2 or greater
  ###  this gets the topographic colors used in other programs by lees
    
    if(missing(ccol)) {   ccol=settopocol()    }
####    program to duplicate (sort of) the the demcmap function in MATLAB
######   n = number of colors

    
    Rz = range(ZTOPO , na.rm = TRUE)

    zee = seq(from=Rz[1], to=Rz[2], length=n)

    a1 = findInterval(zee ,  ccol$calcol$z1)

    cbind(ccol$calcol$z1, ccol$calcol$z2)


    b1 = a1 -1

    cbind(a1, ccol$calcol$z1[a1] , ccol$calcol$z2[a1], zee, (zee-ccol$calcol$z1[a1])/(ccol$calcol$z2[a1]-ccol$calcol$z1[a1])  )


    vecD = (zee-ccol$calcol$z1[a1])/(ccol$calcol$z2[a1]-ccol$calcol$z1[a1])


    R1 = (ccol$calcol$r1[b1]+  vecD*(ccol$calcol$r2[b1]-ccol$calcol$r1[b1]))/255
    G1 = (ccol$calcol$g1[b1]+  vecD*(ccol$calcol$g2[b1]-ccol$calcol$g1[b1]))/255
    B1 = (ccol$calcol$b1[b1]+  vecD*(ccol$calcol$b2[b1]-ccol$calcol$b1[b1]))/255

    R1[R1>1] = 1
    G1[G1>1] = 1
    B1[B1>1] = 1

    tcols = rgb(R1, G1, B1)

    return(tcols)

  }
