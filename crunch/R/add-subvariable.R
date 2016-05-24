##' Add subvariable to an array
##'
##' This function conceals the dirty work in making this happen. The array
##' gets unbound and rebound into a new array with the new variable added.
##' @param variable the array variable to modify
##' @param subvariable the subvariable to add
##' @return a new version of \code{variable} with the indicated subvariables
##' @export
addSubvariable <- function (variable, subvariable){
    ## Store some metadata up front
    payload <- copyVariableReferences(variable)
    subvars <- subvariables(variable)
    subvar.urls <- urls(subvars)
    subvar.names <- names(subvars)

    # TODO: could support taking a VariableDefinition for subvariable
    # if (inherits(subvariable, 'VariableDefinition')) {
    #     ds <- addVariables(ds, subvariable)
    #     subvariable <- ds[[subvariable$alias]]
    # }

    ## Unbind
    old.subvar.urls <- unlist(unbind(variable))

    ## Add the new variable URL to those we had before
    payload$subvariables <- I(c(subvar.urls, self(subvariable)))
    class(payload) <- "VariableDefinition"

    ## Rebind
    new_url <- POSTNewVariable(variableCatalogURL(variable), payload)

    ## Prune subvariable name prefix, or otherwise reset the names
    subvars <- Subvariables(crGET(absoluteURL("subvariables/", new_url)))
    subvar.urls <- c(subvar.urls, self(subvariable))
    subvar.names <- c(subvar.names, name(subvariable))
    names(subvars) <- subvar.names[match(urls(subvars), subvar.urls)]

    ## What to return? This function is kind of a hack.
    invisible(new_url)
}

# addSubvariable <- addSubvariables
