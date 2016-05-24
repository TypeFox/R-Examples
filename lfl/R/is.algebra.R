is.algebra <- function(a) {
    return(is.list(a) && 
           is.function(a$n) &&
           is.function(a$t) &&
           is.function(a$pt) &&
           is.function(a$c) &&
           is.function(a$pc) &&
           is.function(a$r) &&
           is.function(a$b))
}
