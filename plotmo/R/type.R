# type.R: plotmo functions for getting the default type arg for predict() and residuals()
#         this is used when plotmo's argument "type" is NULL (the default)

# get the type for predict()
plotmo.type <- function(object, ...)
{
    UseMethod("plotmo.type")
}
plotmo.type.default <- function(object, ...)
{
    "response"
}
plotmo.type.nnet <- function(object, ...)
{
    "raw"
}
plotmo.type.knn3 <- function(object, ...)
{
    "prob"
}
plotmo.type.tree <- function(object, ...) # tree package
{
    "vector"
}
plotmo.type.fda <- function(object, ...) # mda package
{
    "class"
}
# get the type for residuals()
plotmo.residtype <- function(object, ...)
{
    UseMethod("plotmo.residtype")
}
plotmo.residtype.default <- function(object, ...)
{
    plotmo.type(object, ...) # use the predict type
}
