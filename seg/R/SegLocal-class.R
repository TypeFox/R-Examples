# ------------------------------------------------------------------------------
# Class 'SegLocal'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
setClass(Class = "SegLocal", 
         slots = c(coords = "matrix", data = "matrix", env = "matrix",
                   proj4string = "CRS"))

setValidity(Class = "SegLocal", 
  method = function(object) {
    if (ncol(object@coords) != 2 || !is.numeric(object@coords))
      paste("'coords' must be a numeric matrix with two columns, x and y")
    else if (ncol(object@data) < 2 || !is.numeric(object@data))
      paste("'data' must be a numeric matrix with at least two columns")
    else if (nrow(object@data) != nrow(object@coords))
      paste("'data' must have the same number of rows as 'coords'")
    else if (any(dim(object@env) != dim(object@data)))
      paste("'env' must be a matrix with the same dimensions as 'data'")
    else if (!is.numeric(object@env))
      paste("'env' must be a numeric matrix")
    else if (class(object@proj4string) != "CRS")
      paste("'proj4string' is not a valid CRS object")
    else
      TRUE
  })

SegLocal <- function(coords, data, env, proj4string = CRS(as.character(NA))) {
  new("SegLocal", coords = coords, data = data, env = env,
      proj4string = proj4string)
}

# SegLocal <- function(coords, data, env, proj4string = CRS(as.character(NA))) {
#   # Validate argument 'coords' -----------------------------------------------
#   if (missing(coords))
#     stop("'coords' is missing, with no default", call. = FALSE)
#   else if (!is.matrix(coords)) {
#     coords <- try(as.matrix(coords))
#     if (class(coords) == "try-error")
#       stop("'coords' cannot be coerced to a matrix object", call. = FALSE)
#   }
#   if (ncol(coords) != 2 || !is.numeric(coords))
#     stop("'coords' must be a numeric matrix with two columns", call. = FALSE)
#   
#   # Validate argument 'data' -------------------------------------------------
#   if (missing(data))
#     stop("'data' is missing, with no default", call. = FALSE)
#   else if (!is.matrix(data)) {
#     data <- try(as.matrix(data))
#     if (class(data) == "try-error")
#       stop("'data' cannot be coerced to a matrix object", call. = FALSE)
#   }
#   if (ncol(data) < 2 || !is.numeric(data))
#     stop("'data' must be a numeric matrix with at least two columns", 
#          call. = FALSE)
#   else if (nrow(data) != nrow(coords))
#     stop("'data' must have the same number of rows as 'coords'", 
#          call. = FALSE)
#   
#   # Validate argument 'env' --------------------------------------------------
#   if (missing(env))
#     env <- data
#   else if (!is.matrix(env)) {
#     env <- try(as.matrix(env))
#     if (class(env) == "try-error")
#       stop("'env' cannot be coerced to a matrix object", call. = FALSE)
#   }
#   if (any(dim(env) != dim(data)))
#     stop("'env' must be a matrix with the same dimensions as 'data'", 
#          call. = FALSE)
#   else if (!is.numeric(env))
#     stop("'env' must be a numeric matrix", call. = FALSE)
#   
#   # Validate argument 'proj4string' ------------------------------------------
#   if (class(proj4string) != "CRS")
#     stop("'proj4string' is not a valid CRS object", call. = FALSE)
#   
#   new("SegLocal", coords = coords, data = data, env = env,
#       proj4string = proj4string)
# }
