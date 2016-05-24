roll.rgl.die <- function( side=sample(6,1), steps=250 ) {

  rgl::rgl.viewpoint(45,45)


  tmp <- seq(45, by=90, length=4)
  tmp2 <- c(-1,1,-1,1)
  
  for (j in 1:4) {
    for (i in seq(0,90,length=steps)) {
      rgl::rgl.viewpoint(tmp[j]+i, -tmp2[j]*45+tmp2[j]*i)
    }
  }

  if( side==1 ){
    for(i in seq(0,45, length=steps/2)) {
      rgl::rgl.viewpoint(45+i, 45+i)
    }
  } else if( side==6 ) {
    for(i in seq(0,90, length=steps)) {
      rgl::rgl.viewpoint(45+i, 45-i)
    }
    for(i in seq(0,45, length=steps/2)) {
      rgl::rgl.viewpoint(135+i, -45-i)
    }
  } else {
    tmp3 <- c(NA,3,0,2,1)[side]
    for(j in seq(1,length=tmp3)){
      for(i in seq(0,90,length=steps)) {
        rgl::rgl.viewpoint(tmp[j]+i, -tmp2[j]*45+tmp2[j]*i)
      }
    }
    for(i in seq(0,45, length=steps/2)) {
      rgl::rgl.viewpoint(tmp[tmp3+1]+i, -tmp2[tmp3+1]*45+tmp2[tmp3+1]*i)
    }
  }
  return(side)
}
