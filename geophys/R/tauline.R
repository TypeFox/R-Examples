tauline<-function(Rp, P1, P2, Rview, ES, NN)
  {
    ############ Rp = rotated points describing plane
    ############  P1 and P2 two points extracted from screen
    #####     Rview rotation matrix for viewing
    ######   ES = eigen value decomposition from eigen
    #########   NN normal vector to plan in unrotated coordinates
    ##########   given a 3D perspective plot of a rotated plane
    ###  clicking twice on the figure,
    ###  calculate the shear stress along the line in the plane
    ###   relative to the principle stresses


    ####   first find the normal in the rotated coordinate system:

  RHS = Rp[,3]

    tR = t(Rp)

    A = cbind(t(tR[1:2, 1:3]), rep(1, length=3))

    ATA  = t(A) %*% A
  M = solve(ATA) %*% t(A) %*% RHS



  L = list(x=c(P1$x, P2$x), y=c(P1$y, P2$y))

  
    ZEES = L$x*M[1] + L$y*M[2] + M[3]

    ZPOINTS = cbind(L$x, L$y, ZEES)
    INV = solve(Rview)


    LPts = cbind(ZPOINTS, rep(0, length=length(L$x)))

    LPOINTS = LPts%*% INV

########   un rotation matrix
   
    ###  get vector from one point to the other

    VecLine =  LPOINTS[2 , 1:3] - LPOINTS[1 , 1:3] 
    VecLine = VecLine/sqrt(sum(VecLine*VecLine))

    ########   Dot product:  T %*% SIGMA %*% N

    TAUline =  NN[1]*VecLine[1]*ES$values[1] + NN[2]*VecLine[2] * ES$values[2]  +NN[3]*VecLine[3] *  ES$values[3]

    ###   return the shear stress along this direction
    return(list(tau=TAUline, vec=VecLine))

  }




