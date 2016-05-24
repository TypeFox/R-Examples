`ranef.randomLCA` <-
function(object,...) {
    if (!inherits(object, "randomLCA"))
        stop("Use only with 'randomLCA' objects.\n")
    if (!object$random)
         stop("Object must be fitted with random=TRUE.\n")
    return(object$ranef)
 }

ranef <-
  ## Short form for generic function for extracting the random effects
  function(object, ...) UseMethod("ranef")
