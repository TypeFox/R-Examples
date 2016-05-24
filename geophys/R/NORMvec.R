NORMvec <-
function(PPs, xscale, Rview, aglyph=list(), add=TRUE)
  {
    if(missing(add)) {  add=TRUE   }
    if(missing(aglyph)) {
      headlen =xscale* .3/6
      len =xscale* .7/6
      basethick =xscale* 0.05/2
      headlip =xscale* .02/2
      aglyph = RFOC::Z3Darrow(len = len , basethick =basethick , headlen =headlen , headlip=headlip )
    }

    
    g1  = PPs[1,1:3] - PPs[3, 1:3]
    g2  = PPs[2,1:3] - PPs[3, 1:3]
    g1  = g1/sqrt(sum(g1*g1))
    g2 = g2/sqrt(sum(g2*g2))
    B = unlist(AXB.prod(list(x=g1[1], y=g1[2], z=g1[3]), list(x=g2[1], y=g2[2], z=g2[3])));

    B= B/sqrt(sum(B*B) )
    ### B is now normalized

    L   = list(x1 =PPs[3,1] ,
      y1 =PPs[3,2] ,
      z1=PPs[3,3] ,
      x2 =PPs[3,1]+xscale*B[1]/5 ,
      y2 =PPs[3,2]+xscale*B[2] /5,
      z2=PPs[3,3]+xscale*B[3] /5 
      ) 
    if(add) RFOC::BOXarrows3D(L$x1,L$y1,L$z1, L$x2,L$y2,L$z2,  aglyph=aglyph,  Rview=Rview, col='green')

    return(B)

  }

