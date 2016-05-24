enter <-
function(argument, example = "none", what = "", ...){
    cat("\n\tWhat is the argument ", mark(argument, F),
        " ? - ex: ", mark(example, F), "\n\t\t")
    scan(what = what, ...)
}
