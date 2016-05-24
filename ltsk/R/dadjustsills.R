dadjustsills <-
function(sill0,sills0,sill0t)
{ 
  ## sills : sill(0,0),sill(s,0) and sill(0,t)
  ## values : adjusted sills satisfying p.d. constraints
  d <- sill0 - sills0 - sill0t
  d1 <- sill0 - sills0
  d2 <- sill0 - sill0t
  chk <- any(d1 <0, d2 <0)
  if( chk){  
	  sills0 <- min(sills0,sill0)
	  sill0t <- min(sill0t,sill0)
	  return(c(sill0,sills0,sill0t))
  }
  if (d >=0){
	  sills0 <- sills0 + min(0.51*(d+1e-3),d1)
	  sill0t <- sill0t + min(0.51*(d+1e-3),d2)
	  return(c(sill0,sills0,sill0t))
  }
  c(sill0,sills0,sill0t)
}
