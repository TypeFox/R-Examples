update.sbchoice <- function(object, new, evaluate = TRUE, ...)
{
    call <- getCall(object)

    if (!missing(new)) {
        call$formula <- update(object$formula, new)
    }

    if (evaluate == TRUE) {
        eval(call, parent.frame())
    } else {
        call
    }
}
