if(!isGeneric("covMatrix")) {
  setGeneric(name    = "covMatrix",
             def     = function(object, X, noise.var=NULL) standardGeneric("covMatrix")
  )
}

if(!isGeneric("covMat1Mat2")) {
  setGeneric(name    = "covMat1Mat2",
             def     = function(object, X1, X2, nugget.flag=FALSE) standardGeneric("covMat1Mat2")
  )
}

if(!isGeneric("covparam2vect")) {
  setGeneric(name    = "covparam2vect",
             def     = function(object) standardGeneric("covparam2vect")
             ##  ,package = "DiceKriging"
  )
}

if(!isGeneric("vect2covparam")) {
  setGeneric(name    = "vect2covparam",
             def     = function(object, param) standardGeneric("vect2covparam")
             ##  ,package = "DiceKriging"
  )
}

if(!isGeneric("coef")) {
  setGeneric("coef",
           function(object, ...) standardGeneric("coef")
  )
}

# POUR EVITER LES CONFLITS AVEC gplab
# if(!isGeneric("coef<-")) {
#   setGeneric("coef<-",
#            function(object, ..., value) standardGeneric("coef<-")
#   )
# }

if(!isGeneric("covParametersBounds")) {
  setGeneric(name    = "covParametersBounds",
             def     = function(object, X) standardGeneric("covParametersBounds")
             ##  ,package = "DiceKriging"
  )
}

if(!isGeneric("paramSample")) {
  setGeneric(name    = "paramSample",
             def     = function(object, n, ...) standardGeneric("paramSample")
             ##  ,package = "DiceKriging"
  )
}

if(!isGeneric("covMatrixDerivative")) {
  setGeneric(name = "covMatrixDerivative",
             def = function(object, X, C0, k, ...) standardGeneric("covMatrixDerivative")
  )
}

if(!isGeneric("covVector.dx")) {
  setGeneric(name = "covVector.dx",
             def = function(object, x, X, c) standardGeneric("covVector.dx")
  )
}

if(!isGeneric("show")) {
  setGeneric(name    = "show",
             def     = function(object) standardGeneric("show")
  )
}


