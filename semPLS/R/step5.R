step5 <-
function(Wold, Wnew, tol, converged, convCrit){
  if(convCrit=="relative") diff <- abs((Wold-Wnew)/Wnew) # used
  if(convCrit=="square") diff <- (Wold-Wnew)^2
  
  if (all(diff[!is.nan(diff)] < tol)) converged <- TRUE
  else Wold <- Wnew
  return(list(Wold=Wold, converged=converged))
}
