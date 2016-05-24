##' @rdname makeArray
##' @export
makeMR <- function (list_of_variables, dataset=NULL, pattern=NULL, key=namekey(dataset), name, selections, ...) {
    Call <- match.call(expand.dots=FALSE)

    if (missing(name)) {
        halt("Must provide the name for the new variable")
    }
    if (missing(selections)) {
        halt(paste("Must provide the names of the",
            "category or categories that indicate the dichotomous",
            "selection"))
    }

    Call[[1L]] <- as.name("prepareBindInputs")
    x <- eval.parent(Call)

    ## Get the actual variables so that we can validate
    vars <- lapply(x$variable_urls,
        function (u) CrunchVariable(allVariables(x$dataset)[[u]]))
    are.categorical <- vapply(vars, is.Categorical, logical(1))
    if (!all(are.categorical)) {
        varnames <- vapply(vars[!are.categorical],
            function (x) name(x),
            character(1))
        halt(serialPaste(varnames),
            " are not Categorical variables. Convert them to ",
            "Categorical before combining to Multiple Response")
    }

    ## Validate selections before binding
    catnames <- unique(unlist(lapply(vars,
        function (y) names(categories(y)))))
    if (!all(selections %in% catnames)) {
        halt("Selection(s) not found in variable's categories. ",
            "Category names are: ", serialPaste(catnames))
        ## Could return more useful messaging here
    }

    return(VariableDefinition(subvariables=I(x$variable_urls), name=name,
        type="multiple_response", selected_categories=I(selections), ...))
}
