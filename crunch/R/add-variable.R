##' Add multiple variables to a dataset
##'
##' This function lets you add more than one variable at a time to a dataset.
##' If you have multiple variables to add, this function will be faster than
##' doing \code{ds$var <- value} assignment because it doesn't refresh the
##' dataset's state in between variable POST requests.
##' @param dataset a CrunchDataset
##' @param ... \code{\link{VariableDefinition}}s or a list of
##' VariableDefinitions.
##' @return \code{dataset} with the new variables added (invisibly)
##' @export
##' @importFrom httpcache halt
addVariables <- function (dataset, ...) {
    var_catalog_url <- shojiURL(dataset, "catalogs", "variables")

    ## Get vardefs and validate
    vardefs <- list(...)
    ## Check for whether a list of vardefs passed
    if (length(vardefs) == 1 &&
        is.list(vardefs[[1]]) &&
        !inherits(vardefs[[1]], "VariableDefinition")) {

        vardefs <- vardefs[[1]]
    }
    ## Check that all are VariableDefinitions
    are.vardefs <- vapply(vardefs, inherits, logical(1),
        what="VariableDefinition")
    if (!all(are.vardefs)) {
        halt("Must supply VariableDefinitions")
    }
    ## Check that if values specified, they have the right length
    vardefs <- lapply(vardefs, validateVarDefRows,
        numrows=getNrow(dataset, filtered=FALSE))

    ## Upload one at a time.
    ## TODO: Server should support bulk insert
    new_var_urls <- lapply(vardefs,
        function (x) try(POSTNewVariable(var_catalog_url, x), silent=TRUE))
        ## Be silent so we can throw the errors together at the end

    ## Check for errors
    errs <- vapply(new_var_urls, is.error, logical(1))
    if (any(errs)) {
        if (length(errs) == 1) {
            ## Just one variable added. Throw its error.
            rethrow(new_var_urls[[1]])
        }
        halt("The following variable definition(s) errored on upload: ",
            paste(which(errs), collapse=", "), "\n",
            paste(unlist(lapply(new_var_urls[errs], errorMessage)), sep="\n"))
        ## Could make better error message, return the URLs of the variables
        ## that errored so that user can delete them and start over, etc.
    }

    dataset <- refresh(dataset)
    invisible(dataset)
}

validateVarDefRows <- function (vardef, numrows) {
    ## Pre-check the column length being sent to the server to confirm that
    ## the number of rows matches what's already in the dataset
    if (!any(c("expr", "subvariables") %in% names(vardef))) {
        new <- length(vardef$values)
        old <- numrows ## Just for naming clarity
        if (new == 1 && old > 1) {
            vardef$values <- rep(vardef$values, old)
            new <- old
        }
        if (new == 0) {
            warning("Adding variable with no rows of data", call.=FALSE)
        } else if (old > 0 && new != old) {
            halt("replacement has ", new, " rows, data has ", old)
        }
    }
    return(vardef)
}

POSTNewVariable <- function (catalog_url, variable) {

    do.POST <- function (x) crPOST(catalog_url, body=toJSON(x, digits=15))

    if (!("expr" %in% names(variable))) {
        ## If deriving a variable, skip this and go straight to POSTing
        if (variable$type %in% c("multiple_response", "categorical_array")) {
            ## Assumes: array of subvariables included, and if MR, at least one
            ## category has selected: TRUE, or "selected_categories" given
            if (!("subvariables" %in% names(variable))) {
                halt("Cannot create array variable without specifying",
                    " subvariables")
            }
            ## Three options supported:
            ## (1) Bind together existing subvariables
            ## (2) Create array from single array definition
            ## (3) Create subvariables individually and then bind them
            ## Sniff to see which case we have. If (1) or (2), proceed normally
            is_catvardef <- function (x) {
                all(c("categories", "values") %in% names(x))
            }
            is_binddef <- is.character(variable$subvariables) &&
                !("categories" %in% names(variable))
            is_arraydef <- is_catvardef(variable) &&
                !any(vapply(variable$subvariables, is_catvardef, logical(1)))
            case3 <- !(is_binddef | is_arraydef)
            if (case3) {
                lapply(variable$subvariables, function (x) {
                    Categories(data=x$categories) ## Will error if invalid
                })

                ## Upload subvars, then bind
                var_urls <- lapply(variable$subvariables,
                    function (x) try(do.POST(x)))
                errs <- vapply(var_urls, is.error, logical(1))
                if (any(errs)) {
                    # Delete subvariables that were added, then raise
                    lapply(var_urls[!errs], function (x) crDELETE(x))
                    halt("Subvariables errored on upload")
                }
                # Else prepare to POST array definition
                variable$subvariables <- I(unlist(var_urls))
            }
        }
        if ("categories" %in% names(variable)) {
            Categories(data=variable$categories) ## Will error if cats are invalid
        }
    }
    out <- do.POST(variable)
    invisible(out)
}
