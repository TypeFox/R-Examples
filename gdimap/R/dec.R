
dec <-
function(depth=4, new=TRUE)
{
  s2 <- s2tessel.zorder(depth=depth, viewgrid=FALSE)
  tc.surf <- s2$tcsurf
  pc <- s2$pc
  pc <- pc / max(abs(pc))
  w <- t(apply(abs(pc), 1, norm01))
  ck <- rgb(w)
  pos <- c(0,0,0)
  pcsurf <- cbind(pc[tc.surf,1] + pos[1],
              pc[tc.surf,2] + pos[2] , pc[tc.surf,3])
  # open3d()
  if(new)
   rgl.open()
  rgl.triangles(pc[tc.surf,1], pc[tc.surf,2] ,pc[tc.surf,3],
    col=ck[tc.surf], alpha=0.7)
  # points3d(pc, col=ck)
  ## lines
  z1 <- rbind(c(-1, 0, 0), c(1, 0, 0), c(0,-1, 0), c(0, 1, 0),
    c(0, 0, -1), c(0, 0, 1))
  z2 <- z1 * 1.5
  z <- t(apply(abs(z2), 1, norm01))
  ck <- rgb(z)
  rgl.points(z2, col=ck, alpha=1)
  rgl.lines(z2, col=ck, lwd=2, alpha=1)
  rgl.texts(c(2,0,0), text="X", col="red")
  rgl.texts(c(0,2,0), text="Y", col="green")
  rgl.texts(c(0,0,2), text="Z", col="blue")
  # rgl.viewpoint(0,0)
  rgl.viewpoint(30,30)
  par3d("zoom"=0.6)
}


