DBridge <- function(x=0, y=0, t0=0, T=1, delta, drift, sigma, ...){
 done <- FALSE
 while(!done){
  Y1 <- sde.sim(X0=x, drift=drift, sigma=sigma, t0=t0, T=T, delta=delta, ...)
  Y2 <- sde.sim(X0=y, drift=drift, sigma=sigma, t0=t0, T=T, delta=delta, ...)
  Y3 <- ts(rev(Y2), start=start(Y2), end=end(Y2),deltat=deltat(Y2))

  id <- Inf
  if(Y1[1]>=Y3[1]){
   if(!all(Y1>Y3))
    min(which(Y1 <= Y3))-1 -> id
   } else {
    if(!all(Y1<Y3))
     min(which(Y1 >= Y3))-1 -> id
   }

  if(id==0 || id==length(Y1) || id == Inf) {
    message("\nno crossing, trying again...")
	done <- FALSE
  } else {
   done <- TRUE
  }
 }
 
  B <- ts(c(Y1[1:id], Y3[-(1:id)]), start=start(Y1),end=end(Y1),frequency=frequency(Y1))
  return(invisible(B))
}
