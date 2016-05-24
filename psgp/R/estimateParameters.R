estimateParameters.psgp = function(object, ...) {

  origObs = object$observations
  
  rotated = FALSE
  if (object$params$doAnisotropy) {
      object = estimateAnisotropy(object) 
      #rotate Data
    if (object$anisPar$doRotation && all(as.character(object$formulaString[[3]])=="1"))
      object$observations=rotateAnisotropicData(object$observations,object$anisPar)
    rotated = TRUE
  }  

  
  # Restore measurement error characteristics
  # object$observations$oeid = origObs$oeid;
  # object$observations$oevar = origObs$oevar
  
  # Estimate parameters using PSGP
  object = learnParameters(object)

  # Restore original observations
  if (rotated) 
    object$observations = origObs
  object
}