#
# whole sphere tesselation  
#
s2tessel.zorder <-
function(depth=3, viewgrid=FALSE, saveg=FALSE)
{
  l1 <- icosahedron3d(col="magenta")
	l2 <- subdivision3d( l1, depth=depth)
	g0 <- t(l2$vb)
	g0 <- g0[,1:3]/g0[,4]
	g0 <- as.matrix(g0)
	g0 <- g0/max(g0)
	## sort in ZZ and use z-up-and-down
	zi <- sort(g0[,3], ind=TRUE, decreasing=TRUE)
	g1 <- g0[zi$ix,]
	nz1 <- which(g1[,3] > 0)
	nz2 <- which(g1[,3] <= 0)
	g0 <- rbind(g1[nz1,], g1[nz2,])
	if(saveg) {
		f <- tempfile(pattern="g0",fileext=".Rdata")
		write.table(g0, file=f, row.names=FALSE, col.names=FALSE)
		cat("wrote",f,"\n")
	}
	tc <-  geometry::delaunayn(g0)
	tc.surf <- t( surf.tri(g0,tc) )
	s2tess <- list(pc=g0, tcsurf=tc.surf)
	#---------
	if(viewgrid) {
		open3d()
		plot3d(g0) 
		rgl.triangles(g0[tc.surf,1], g0[tc.surf,2] ,g0[tc.surf,3], col="blue", alpha=0.2)
 		for(i in 1:(8*360)) rgl.viewpoint(i/8)
	}
	invisible(s2tess)
}

