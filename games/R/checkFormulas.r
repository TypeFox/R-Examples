##
## INPUT:
## f: object inheriting from class "formula", or a list of such objects
## argname: character string specifying the name of the argument being checked
## in the original function (in order to give an informative error message in
## case of failure)
##
## RETURN:
## object of class "Formula", combining supplied formulas (if 'f' is a list)
## into a big one with multiple right-hand sides
##
checkFormulas <- function(f, argname = "formulas")
{
    if (inherits(f, "list")) {
        f <- do.call(as.Formula, unname(f))
    } else if (inherits(f, "formula")) {
        f <- as.Formula(f)
    } else {
        stop(argname, " must be a list of formulas or a formula")
    }

    return(f)
}
