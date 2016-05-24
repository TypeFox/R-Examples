rgl.die <- function(x = 1:6, col.cube='white',col.pip='black',sides=x, ...) {

  if(!requireNamespace('rgl', quietly = TRUE)) stop("This function depends on the 'rgl' package wich is not available")

  rgl::rgl.viewpoint(45,45)

  pip.coords <- function( x,y ) {
    xc <- yc <- numeric(0)
    for(i in 0:39){
      xc <- c(xc, x, 0.05*cos(pi/20*i)+x, 0.05*cos(pi/20*(i+1))+x)
      yc <- c(yc, y, 0.05*sin(pi/20*i)+y, 0.05*sin(pi/20*(i+1))+y)
    }
    cbind(xc,yc)
  }

  pip.loc <- list(matrix( 0.5, ncol=2, nrow=1),
                  cbind( c(.25, .75), c(.25, .75)),
                  cbind( c(.25, .5, .75), c(.25, .5, .75)),
                  cbind( c(.25, .25, .75, .75), c(.25, .75, .75, .25)),
                  cbind( c(.25, .25, .75, .75, .5), c(.25, .75, .75, .25, .5)),
                  cbind( c(.25, .25, .25, .75, .75, .75),
                         c(.25, .5, .75, .75, .5, .25)))

  rgl::rgl.quads( c(0,0,1,1), c(0,1,1,0), c(0,0,0,0), col=col.cube)
  rgl::rgl.quads( c(0,0,1,1), c(0,1,1,0), c(1,1,1,1), col=col.cube)
  rgl::rgl.quads( c(0,0,0,0), c(0,1,1,0), c(0,0,1,1), col=col.cube)
  rgl::rgl.quads( c(1,1,1,1), c(0,1,1,0), c(0,0,1,1), col=col.cube)
  rgl::rgl.quads( c(0,0,1,1), c(0,0,0,0), c(0,1,1,0), col=col.cube)
  rgl::rgl.quads( c(0,0,1,1), c(1,1,1,1), c(0,1,1,0), col=col.cube)

  tmp <- pip.loc[[ sides[1] ]]
  for( i in 1:nrow(tmp) ){
    xy <- pip.coords( tmp[i,1], tmp[i,2] )
    rgl::rgl.triangles(xy[,1], rep(1.001, nrow(xy)), xy[,2], col=col.pip,lit=FALSE)
  }

  tmp <- pip.loc[[ sides[2] ]]
  for( i in 1:nrow(tmp) ){
    xy <- pip.coords( tmp[i,1], tmp[i,2] )
    rgl::rgl.triangles(xy[,1], xy[,2], rep(1.001, nrow(xy)), col=col.pip,lit=FALSE)
  }

  tmp <- pip.loc[[ sides[3] ]]
  for( i in 1:nrow(tmp) ){
    xy <- pip.coords( tmp[i,1], tmp[i,2] )
    rgl::rgl.triangles( rep(1.001, nrow(xy)), xy[,1], xy[,2], col=col.pip,lit=FALSE)
  }

  tmp <- pip.loc[[ sides[4] ]]
  for( i in 1:nrow(tmp) ){
    xy <- pip.coords( tmp[i,1], tmp[i,2] )
    rgl::rgl.triangles( rep(-0.001, nrow(xy)), xy[,1], xy[,2], col=col.pip,lit=FALSE)
  }

  tmp <- pip.loc[[ sides[5] ]]
  for( i in 1:nrow(tmp) ){
    xy <- pip.coords( tmp[i,1], tmp[i,2] )
    rgl::rgl.triangles(xy[,1], xy[,2], rep(-0.001, nrow(xy)), col=col.pip,lit=FALSE)
  }

  tmp <- pip.loc[[ sides[6] ]]
  for( i in 1:nrow(tmp) ){
    xy <- pip.coords( tmp[i,1], tmp[i,2] )
    rgl::rgl.triangles(xy[,1], rep(-0.001, nrow(xy)), xy[,2], col=col.pip,lit=FALSE)
  }

}
