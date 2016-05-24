TriangleCenter<-function(P1, P2, P3, A1= 0, A2= 360, KNum=10 )
  {

    ###  given three points in 3D space find the
      ###    center of the triangle formed by these

    if(missing(A1)) { A1 = 0 }
    if(missing(A2)) { A2 = 360 }
    if(missing(KNum)) { KNum=10 }

    

    ampvec<-function(a){return(sqrt(sum(a*a))) }

    check = FALSE
    
    
    if(is.list(P1)) P1 = unlist(P1)
    if(is.list(P2)) P2 = unlist(P2)
    if(is.list(P3)) P3 = unlist(P3)

    ###  conbine the points into a matrix
    PP = rbind(P1, P2, P3, rep(1, 3) )

    ####   get vector from point 1 to 3 and from point 2 to 3
    g13 = P1-P3
    g23 = P2-P3

    ####  get normal to the 3D triangle
    B = unlist(AXB.prod(g13, g23))
    ###  normalize
    B = B/sqrt(sum(B * B))

####  find the rotation that will translate the triangle
 ####     from it position in space to the the origin where
 ####     it lies in the XY-plane
####  the point 3 is the anchor - or origin of hte plot in 2D
    
MF = rot2Zplane(B, P3) 

PP = cbind( rbind(P1, P2, P3),  rep(1, 3))

    

  
UMAT = (PP) %*% MF

######  check here ---- problems

    if(check)
      {
        dev.set(dev.next())
        plot(UMAT[,1], UMAT[,2], asp=1)
        text(UMAT[,1], UMAT[,2], labels=1:3, pos=1, xpd=TRUE)
        
        u13 =  c( (UMAT[1,1]-UMAT[3,1]),  (UMAT[1,2]-UMAT[3,2]))
        u23 =  c( (UMAT[2,1]-UMAT[3,1]),  (UMAT[2,2]-UMAT[3,2]))

        dotu12  = sum( u13*u23)/(ampvec(u13)*ampvec(u23))
        180*acos(dotu12)/pi

        
        dotg12 = sum(g13 * g23)/(ampvec(g13)*ampvec(g23))
        180*acos(dotg12)/pi

        arrows(UMAT[3,1], UMAT[3,2], UMAT[3,1]+u13[1], UMAT[3,2]+u13[2])
        arrows(UMAT[3,1], UMAT[3,2], UMAT[3,1]+u23[1], UMAT[3,2]+u23[2])

        
      }

    AA = TriangleInfo(UMAT[1,1:2], UMAT[2,1:2], UMAT[3,1:2], add=check)   
    if(check)   dev.set(dev.prev())
    
cc = c( AA$BI, 0.0, 1)


flatcirc =   GEOmap::darc(rad = AA$r, ang1 = A1, ang2 = A2, x1 = cc[1], y1 = cc[2], n = KNum) 

cr = cc %*% solve(MF)

JC = cbind(flatcirc$x,flatcirc$y, rep(0, length(flatcirc$x)),  rep(1, length(flatcirc$x))  )
 
Cinscribed = JC  %*% solve(MF)


cR = c( AA$CIRCUM, 0.0, 1) %*% solve(MF)
    
    
##   UMAT %*% solve(MF)

  ##  print(cr)
    
    invisible(list(Center=cr, r=AA$r   , Cinscribed=Cinscribed, CIRCUM=cR ) )

    
  }
