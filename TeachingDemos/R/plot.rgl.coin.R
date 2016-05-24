rgl.coin <- function(x, col='black', heads=x[[1]],
                          tails=x[[2]], ... ) {
    if(!requireNamespace('rgl', quietly = TRUE)) stop("This function depends on the 'rgl' library which is not available")

  if(missing(x)) x <- TeachingDemos::coin.faces

  rgl::rgl.viewpoint(0,0)

  for(i in 0:39) {
    rgl::rgl.triangles(c(.5, cos(pi/20*i)/2+0.5, cos(pi/20*(i+1))/2+0.5),
                  c(.5, sin(pi/20*i)/2+0.5, sin(pi/20*(i+1))/2+0.5),
                  c(0,0,0))
    rgl::rgl.triangles(c(.5, cos(pi/20*i)/2+0.5, cos(pi/20*(i+1))/2+0.5),
                  c(.5, sin(pi/20*i)/2+0.5, sin(pi/20*(i+1))/2+0.5),
                  c(0.03,0.03,0.03))
    rgl::rgl.quads( c(cos(pi/20*i)/2+0.5, cos(pi/20*i)/2+0.5,
                 cos(pi/20*(i+1))/2+0.5, cos(pi/20*(i+1))/2+0.5),
               c(sin(pi/20*i)/2+0.5, sin(pi/20*i)/2+0.5,
                 sin(pi/20*(i+1))/2+0.5, sin(pi/20*(i+1))/2+0.5),
               c(0,0.03,0.03,0)
              )
  }

  tmp <- rep( 1:nrow(heads), each=2 )
  tmp <- c(tmp[-1],1)

  rgl::rgl.lines( heads[tmp,1], heads[tmp,2], rep(0.035, length(tmp) ),
             col=col, lit=FALSE)

  tmp <- rep( 1:nrow(tails), each=2 )
  tmp <- c(tmp[-1],1)

  rgl::rgl.lines( tails[tmp,1], tails[tmp,2], rep(-0.005, length(tmp) ),
             col=col, lit=FALSE)

}


#coin.faces <- list( qh=cbind( c(.5,.5), c(.75,.25) ),
#                    qt=cbind( c(.5, .25, .5, .75, .5),
#                      c(.75, .5, .25, .5, .75)) )
