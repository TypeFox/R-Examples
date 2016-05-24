localSCovarianceGradients <-
function(model) {
  gK = -model$d*model$invK_uu + model$invK_uu%*%model$S%*%model$invK_uu
  gK = gK * .5
  return (gK)
}
