se.fixef <- function(model) {
    #stopifnot(is(model, "mer"))
    sqrt(diag(vcov(model)))
}

