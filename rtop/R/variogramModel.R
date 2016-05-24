rtopVariogramModel = function(model = "Ex1", sill = NULL, range = NULL, exp = NULL,
    nugget = NULL, exp0 = NULL,
    observations = NULL, formulaString = obs~1) {
if (model == "ex1") {
  if (!is.null(observations)) {
    parInit = findParInit(formulaString, observations, model)$par0
    if (is.null(sill)) sill = parInit[1]
    if (is.null(range)) range = parInit[2]
    if (is.null(nugget)) nugget = parInit[3]
    if (is.null(exp)) exp = parInit[4]
    if (is.null(exp0)) exp0 = parInit[5]
  } else {
    if (is.null(sill)) sill = 1
    if (is.null(range)) range = 1
    if (is.null(exp)) exp = 0
    if (is.null(nugget)) nugget = 0
    if (is.null(exp0)) exp0 = 1
  }
  variogramModel = list(model = model, params = c(sill, range, nugget, exp, exp0))
  class(variogramModel) = "rtopVariogramModel"
  }
  variogramModel
}

updateRtopVariogram.rtop = function(object, ...) {
object$variogramModel = updateRtopVariogram(object$variogramModel, sampleVariogram = object$variogram,
observations = object$observations, ...)
object
}

updateRtopVariogram.rtopVariogramModel = function(object, action = "mult", ..., checkVario = FALSE, 
sampleVariogram = NULL, observations = NULL){
  variogramModel = object
  dots = list(...)

  if (variogramModel$model == "Ex1") {
    if ("sill" %in% names(dots)) {
      if (action == "mult") {
        variogramModel$params[1] = variogramModel$params[1]*dots$sill
      } else if (action == "add") {
        variogramModel$params[1] = variogramModel$params[1]*dots$sill
      } else if (action == "replace") {
        variogramModel$params[1] = dots$sill
      }
    }
    if ("range" %in% names(dots)) {
      if (action == "mult") {
        variogramModel$params[2] = variogramModel$params[2]*dots$range
      } else if (action == "add") {
        variogramModel$params[2] = variogramModel$params[2]*dots$range
      } else if (action == "replace") {
        variogramModel$params[2] = dots$range
      }
    }
    if ("nugget" %in% names(dots)) {
      if (action == "mult") {
        variogramModel$params[3] = variogramModel$params[3]*dots$nugget
      } else if (action == "add") {
        variogramModel$params[3] = variogramModel$params[3]*dots$nugget
      } else if (action == "replace") {
        variogramModel$params[3] = dots$nugget
      }
    }
    if ("exp" %in% names(dots)) {
      if (action == "mult") {
        variogramModel$params[4] = variogramModel$params[4]*dots$exp
      } else if (action == "add") {
        variogramModel$params[4] = variogramModel$params[4]*dots$exp
      } else if (action == "replace") {
        variogramModel$params[4] = dots$exp
      }
    }
    if ("exp0" %in% names(dots)) {
      if (action == "mult") {
        variogramModel$params[5] = variogramModel$params[5]*dots$exp0
      } else if (action == "add") {
        variogramModel$params[5] = variogramModel$params[5]*dots$exp0
      } else if (action == "replace") {
        variogramModel$params[5] = dots$exp0
      }
    }
    if (checkVario) checkVario(variogramModel, sampleVariogram = sampleVariogram, observations = observations)  
  }
  variogramModel    
}





plot.rtopVariogramCloud = function(x,  ...) {
x$np = x$ord
class(x) = "variogramCloud"
plot(x, ...)
}
