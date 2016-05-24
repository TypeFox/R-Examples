covparms <-
function(object) {
  if(class(object) != "glmssn") return("Not a glmssn object")
  summary(object)$covariance.parameter.estimates
}

