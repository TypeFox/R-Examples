vcov.dbchoice <- function(object, ...)
{
    solve(object$Hessian)
}
