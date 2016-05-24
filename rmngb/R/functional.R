composition <- function(...) {
    listFunctions <- list(...)
    stopifnot(length(listFunctions) > 1,
              all(vapply(listFunctions, is.function, logical(1))))
    
    comp <- function(l) {
        if (length(l) == 1) {
            l[[1]]
        } else {
            function(...) l[[1]](comp(l[- 1])(...))
        }
    }
    
    return(comp(listFunctions))
}

cons <- function(a, b)
    function(f) f(a, b)
car <- function(x)
    x(function(a, b) a)
cdr <- function(x)
    x(function(a, b) b)
cadr <- function(x)
    car(cdr(x))