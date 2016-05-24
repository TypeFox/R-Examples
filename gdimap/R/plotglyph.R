plotglyph <-
function(odf, grad, pk, kdir=6, vmfglyph=TRUE, pos=c(0,0,0))
{
  sz <- length(odf)
  pc0 <- grad
  # pc0 <- rbind(grad, -grad)
  tc <-  geometry::delaunayn(pc0)
  tc.surf <- t( surf.tri(pc0,tc) )
  #  open3d()
  odf <- odf - min(odf)
  gk <- genfa(odf)
  ##-------------
  ## RGB channels
  # calibrated colors for each odf-vertice
  zch <- pc0 * gk
  zch <- t(apply(abs(zch), 1, norm01)) 
  ck <- rgb(zch)
  ss <- sort(odf, decreasing=TRUE, index=TRUE)
  pc <- pc0 * as.vector(odf) 
  pc <- pc / max(pc)
  # 1st direction of max odf values for odfs
  ix <- ss$ix[1] 
  pcsurf <- cbind(pc[tc.surf,1]+pos[1],
              pc[tc.surf,2]+pos[2] , pc[tc.surf,3])+pos[3]
  # open3d()
  # plot3d(pc, col=ck, xlim=c(-1,1), ylim=c(-1,1), zlim=c(-1,1)) # show color of points
  # points3d(pc)
  # build color for each triangle based on vertices' color 
  rgl.triangles(pc[tc.surf,1]+pos[1], pc[tc.surf,2]+pos[2],
    pc[tc.surf,3]+pos[3], col=ck[tc.surf], alpha=0.3)
  #--------------------
  # visualize actual fiber directions by max(odf)
  np <- pk$np
  if(!np) {
    warning("zero directions")
  }
  na <- min(np,kdir)
  k <-seq(2,2*na, by=2)
  v <- matrix(0, nrow=2*na, ncol=3)
  vpos <- matrix(rep(pos, 2*na), byrow=TRUE, ncol=3)
  if(vmfglyph) {
    pcoords <- pk$pcoords
    v[(1:na)*2, ] <- t(pcoords[,1:na])
  }
  else { # for fpeak glyphs
    v[(1:na)*2-1, ] <- -t(pk$pcoords[,1:na]) 
    v[(1:na)*2, ] <- t(pk$pcoords[,1:na])
  }
  segments3d(v+vpos, add=TRUE, col=rep(1:kdir, each=2), lwd=3, alpha=1)
  rgl.bg(color="white")
  rgl.viewpoint(0,0)
}

