flip.rgl.coin <- function(side=sample(2,1), steps=150) {

  for (i in seq(0,(5+side)*180, length=steps*(5+side)) ){
    rgl::rgl.viewpoint(i,0)
  }
  return(side)
}
