# Methods for lsm.list objects

# Summary method for an lsm.list
summary.lsm.list <- function(object, ...)
    lapply(object, function(x) {
        if (inherits(x, "summary.ref.grid"))  x
        else summary(x, ...)
    })

print.lsm.list <- function(x, ...) 
    print(summary(x, ...))

str.lsm.list = function(object, ...) {
    for(nm in names(object)) {
        cat(paste("$", nm, "\n", sep=""))
        str(object[[nm]])
        cat("\n")
    }
}

# Courtesy methods to make it more friendly for follow-ups
contrast.lsm.list = function(object, ... , which = 1) {
    contrast(object[[which]], ...)
}

pairs.lsm.list = function(x, ..., which = 1) {
    pairs(x[[which]], ...)
}

test.lsm.list = function(object, ..., which = 1) {
    test(object[[which]], ...)
}

confint.lsm.list = function(object, ..., which = 1) {
    confint(object[[which]], ...)
}

cld.lsm.list = function(object, ..., which = 1) {
    cld(object[[which]], ...)
}


