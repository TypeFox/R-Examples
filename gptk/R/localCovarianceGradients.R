localCovarianceGradients <-
function(model, y, dimension) {
  if (!"isSpherical" %in% names(model) || model$isSpherical) {
    invKy = model$invK_uu %*% y
    gK = -model$invK_uu + invKy%*%t(invKy)
  } else {
    if (model$isMissingData)
      m = y[model$indexPresent[[dimension]]]
    else
      m = y

    invKy = model$invK_uu[[dimension]] %*% m
    gK = -model$invK_uu[[dimension]] + invKy%*%t(invKy)
  }
  gK = gK * .5
  return (gK)
}
