does_not_give_warning <- function (regexp = NULL) {
    function(expr) {
        res <- evaluate::evaluate(substitute(expr), parent.frame())
        warnings <- sapply(Filter(evaluate::is.warning, res), "[[", "message")
        if (!is.null(regexp)) {
            matches(regexp, all = FALSE)(warnings)
        }
        else {
            expectation(length(warnings) == 0, "warnings given")
        }
    }
}
