migpdCoefs <-
  # Coerce coefficients of a migpd object.
  # Useful when coefficients are coming from
  # a model with covariates and you want to
  # learn something about the dependence between
  # margins.
function(object, which, coefs){
  if (class(object) != "migpd"){
    stop("object must be of class \'migpd\'")
  }

  if (length(which) != length(coefs)){
    stop("which and coefs should have the same length")
  }

  if (length(which) == 1){
    object$models[[which]]$coefficients <- coefs
  }

  for (i in 1:length(which)){
    object$models[[i]]$coefficients <- coefs[[i]]
  }
  object
}

