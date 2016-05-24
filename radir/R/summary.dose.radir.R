summary.dose.radir <-
function(object, ...)
{
  # Calibrative density statistics
  if (class(object)!="dose.radir") stop("Wrong object")
  mod_cd <- object[[2]][which.max(object[[1]])]
  cd  <- approxfun(object[[2]], object[[1]])
  cd2 <- approxfun(object[[2]], object[[1]]*object[[2]])
  expect_cd <- integrate(cd2, lower=min(object[[2]]), upper=max(object[[2]]))$value
  cd3 <- approxfun(object[[2]], ((object[[2]]-expect_cd)^2)*object[[1]])
  var_cd <- integrate(cd3, lower=min(object[[2]]), upper=max(object[[2]]))$value 
  
  ci <- ci.dose.radir(object, .95)
  
  ans <- data.frame(round(mod_cd,3), round(expect_cd,3), round(sqrt(var_cd),3), paste0("(", round(ci[1],3), "; ", round(ci[2],3), ")"))
  class(ans)  <- "summary.dose.radir"
  return(ans)
}
