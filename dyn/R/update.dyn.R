update.dyn <- function(object, ..., evaluate = TRUE) {
    .Class <- class(object) <- setdiff(class(object), "dyn")
    call <- NextMethod("update", object, ..., evaluate = FALSE)
    class(call$formula) <- c("dyn", class(formula(object)))
    if (evaluate) 
        eval(call, parent.frame())
    else call
}
