# Copied directly from Yplant code (translated from Delphi). 
# Courtesy of Bob Pearcy
madeleafdirection <- function(mor,lor,lar){

  if(lar == 0){  # then  //This condition indicates that the leaf is flat
    ld.e <- sin(mor)
    ld.n <- cos(mor)
    ld.z <- 0
    wd.e <- cos(mor)
    wd.n <- -sin(mor)
    wd.z <- 0
  }  else {
 
    #{get the leaf plane normal vector direction}
    normal.e <- sin(lar)*sin(lor)
    normal.n <- sin(lar)*cos(lor)
    normal.z <- cos(lar)
    #{get the midrib vector direction--leaf plane Y axis}
    if(normal.z == 0){
      ld.e <- 0
      ld.n <- 0
      ld.z <- -1
    } else {
      ld.e <- sin(mor)
      ld.n <- cos(mor)
      ld.z <- -(ld.e*normal.e+ld.n*normal.n)/normal.z
	}
	
    model <- sqrt(ld.e^2 + ld.n^2 + ld.z^2)
    ld.e <- ld.e/model 
	  ld.n <- ld.n/model
	  ld.z <- ld.z/model
    
    #{get leaf plane X axis} -- I think he means leaf width vector, perpendicular to leaf length vector
    wd.e <- ld.n*normal.z - ld.z*normal.n
    wd.n <- ld.z*normal.e - ld.e*normal.z
    wd.z <- ld.e*normal.n - ld.n*normal.e
	}
return(list(ld=c(ld.e,ld.n,ld.z), wd=c(wd.e,wd.n,wd.z)))
}

