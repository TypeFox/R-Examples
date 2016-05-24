profileVarEst <-
function(se0, sg0, D1, D2, r, n1, n2, control = NULL) {
  if (is.null(control)) {
    control <- list(lklMethod = "REML", s2profile = "se", lb = 0.01, ub = 100, tol = 0.01)
  } 
  
  if (control$s2profile == "se"){
      x0 = sg0 / se0 
      S = profile.newton.se(x0, D1, D2, r, n1, n2, control) 
  } else { 
      x0 = se0 / sg0
      S = profile.newton.sg(x0, D1, D2, r, n1, n2, control)
  }
  
  return(list(s2e = S$s2e, s2g = S$s2g, tau = S$tau))
}
