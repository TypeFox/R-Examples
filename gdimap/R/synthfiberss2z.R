synthfiberss2z <-
function(g0, angles=c(20,100), b=3000, S0=1, sigma=NULL, logplot=TRUE, pos=c(0,0,0), showglyph=FALSE, new=TRUE, wi=NULL)
{
  egv <- c(1700, 200, 200) * 10^(-6) # see alexander-2002
  Du <- matrix(egv, ncol=1)
  ## print(Du)
  sv <- numeric(dim(g0)[1])
  na <- length(angles)
  # fi <- 1/na # equal weight for fibers, no isotropic part
	if(is.null(wi)) {
	  l <- length(angles)
	  wi <- rep(1/l, times=l)
	}
	else 
		stopifnot(length(wi) == length(angles))
  for(k in 1:na) {
    angl <- angles[k]
    if(angl) {
      phir <- angl*pi/180 
      grad <- rotate3d(g0, phir, 0, 0, 1)
    }
    else { grad <- g0 }
    gx2 <- grad^2  
    gxx2 <- gx2[,1]; gyy2 <- gx2[,2]; gzz2 <- gx2[,3];
    X <- -b * cbind(gxx2, gyy2, gzz2)
    X <- as.matrix(X)
    y <- X%*%Du
    # y <- fi * exp(y)
    y <- wi[k] * exp(y)
    ## plot(y, pch=20, ty="p")
    ## open3d()
    ## plot3d(ellipse3d(D), col="green", box=FALSE, alpha=0.5)
    sv <- sv + y
  }
  svnfree <- as.vector(sv * S0)
  ## Riccian type noise model of a typical standard deviation leve
  if(!is.null(sigma)) {
    dm <- dim(g0)[1]
     sv <- sqrt((svnfree + sigma*rnorm(dm))^2+(sigma*rnorm(dm))^2)
  }
  else { 
    sv <- svnfree
  }
  if(showglyph) {
    svp <- sv
    ## cat("signal range:", range(sv), "\n")
    if(logplot) 
      svp <- log(svp) # apply to view orientation
    tc <- geometry::delaunayn(g0)
    pc <- g0*svp
    if(new)
      open3d()
    pc <- pc/max(pc)
    tc.surf <- t( surf.tri(pc,tc) )
    lim=c(-1,1); 
    ## plot3d(pc[,1], pc[,2], pc[,3], xlim=lim, ylim=lim, zlim=lim)
    rgl.triangles(
      pc[tc.surf,1]+pos[1],
      pc[tc.surf,2]+pos[2],
      pc[tc.surf,3]+pos[3],
      col="blue", alpha=0.3,
      xlim=lim, ylim=lim, zlim=lim)
    ## visualize fibers 
    v <- matrix(0, nrow=2, ncol=3)
    for(k in 1:na) {
      angl <- angles[k]
      fi <- c(cos(angl*pi/180), sin(angl*pi/180), 0)
      v[1,] <- -fi + pos
      v[2,] <- fi + pos
      segments3d(v, add=TRUE, col="red", lwd=3, alpha=1)
    }
    rgl.viewpoint(0,0)
  }
  invisible(sv)
}

