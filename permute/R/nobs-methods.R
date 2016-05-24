## nobs generic is defined only in R 2.13.0: define here for older R
if (getRversion() < "2.13.0")
    nobs <- function(object, ...) UseMethod("nobs")

## add some nobs() methods - need to be documented
`nobs.numeric` <- function(object, ...) {
    length(object)
}

`nobs.integer` <- function(object, ...) {
    nobs.numeric(object, ...)
}

`nobs.matrix` <- function(object, ...) {
    NROW(object)
}

`nobs.data.frame` <- function(object, ...) {
    NROW(object)
}

`nobs.factor` <- function(object, ...) {
    length(object)
}
