#' @export
funvec <-
function (x, data = parent.frame(), fun = mean, groups = NULL, 
    subset = TRUE, ...) 
{
    dots <- list(...)
    groups <- eval(substitute(groups), data, parent.frame())
    subset <- eval(substitute(subset), data, parent.frame())
    form <- latticeParseFormula(x, data, subset = subset, groups = groups, 
        )
    gm <- tapply(form$left, form$right, fun)
    return(gm[form$right])
}
