
ctm <- function(response, interacting = NULL, shifting = NULL, data = NULL,
                todistr = c("Normal", "Logistic", "MinExtrVal"),
                sumconstr = inherits(interacting, c("formula", "formula_basis")), ...) {

    ### mkgrid() will not work if data is missing
    if (.is.formula(response)) 
        response <- as.basis(response, data = data)
    if (.is.formula(interacting)) 
        interacting <- as.basis(interacting, data = data)
    if (.is.formula(shifting)) 
        shifting <- as.basis(shifting, data = data, remove_intercept = TRUE, ...)

    if (is.character(todistr))
        todistr <- .distr(todistr)

    bases <- list(response = response, interacting = interacting, shifting = shifting)

    if (!is.null(interacting))
        interacting <- b(iresponse = response, iinteracting = interacting, 
                         sumconstr = sumconstr)

    if (is.null(interacting) && is.null(shifting)) {
        mod <- c(bresponse = response)
    } else if (!is.null(interacting) && is.null(shifting)) {
        mod <- c(binteracting = interacting)
    } else if (is.null(interacting) && !is.null(shifting)) {
        mod <- c(bresponse = response, bshifting = shifting)
    } else {
        mod <- c(binteracting = interacting, bshifting = shifting)
    }
    ret <- list(model = mod, response = variable.names(response), 
                todistr = todistr, bases = bases)
    class(ret) <- "ctm"
    return(ret)
}

model.matrix.ctm <- function(object, data, ...)
    return(model.matrix(object$model, data = data, ...))

variable.names.ctm <- function(object, ...)
    variable.names(object$model)

