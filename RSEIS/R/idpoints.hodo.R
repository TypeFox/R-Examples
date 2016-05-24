`idpoints.hodo` <-
function(nbaz, sx, X, Y)
  {

    if(length(X)<1) return(NULL)

    V = nbaz[,1]
    N = nbaz[,2]
    E = nbaz[,3]

    pts = vector()
    
    for(i in 1:length(X))
      {
        plt1 = 1+floor(X[i])
        if(plt1==1)
          {
            x  = RPMG::RESCALE(E, 0, 1, sx[1], sx[2])
            y  = RPMG::RESCALE(N, 0, 1, sx[1], sx[2])
          }

        if(plt1==2)
          {
            x  = RPMG::RESCALE(E, 1, 2, sx[1], sx[2])
            y  = RPMG::RESCALE(V, 0, 1, sx[1], sx[2])
          }

        if(plt1==3)
          {
            x  = RPMG::RESCALE(N, 2, 3 , sx[1], sx[2])
            y  = RPMG::RESCALE(V, 0, 1, sx[1], sx[2])
          }

        dis = (x-X[i])^2+(y-Y[i])^2

        pts[i] = which.min(dis)


      }

    return(pts)
    
  }

