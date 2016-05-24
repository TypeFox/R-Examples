##
##  h y p e r b o l i c . R  Hyperbolic Functions
##


### More trigonometric functions ###

# cotangens
cot <- function(z) {
    stopifnot(is.numeric(z) || is.complex(z))
    1 / tan(z)
}

# cosecans
csc <- function(z) {
    stopifnot(is.numeric(z) || is.complex(z))
    1 / sin(z)
}

# secans
sec <- function(z) {
    stopifnot(is.numeric(z) || is.complex(z))
    1 / cos(z)
}

# arcus cotangens
acot <- function(z) {
    stopifnot(is.numeric(z) || is.complex(z))
    atan(1/z)
}

#arcus cosecans
acsc <- function(z) {
    stopifnot(is.numeric(z) || is.complex(z))
    asin(1/z)
}

# arcus secans
asec <- function(z) {
    stopifnot(is.numeric(z) || is.complex(z))
    acos(1/z)
}


### More hyperbolic functions ###

# hyperbolic cotangens
coth <- function(z) {
    stopifnot(is.numeric(z) || is.complex(z))
    1 / tanh(z)
}

# hyperbolic cosecans
csch <- function(z) {
    stopifnot(is.numeric(z) || is.complex(z))
    1 / sinh(z)  # 2 / (exp(z) - exp(-z))
}

# hyperbolic secans
sech <- function(z) {
    stopifnot(is.numeric(z) || is.complex(z))
    1 / cosh(z)  # 2 / (exp(z) + exp(-z))
}

# area cotangens
acoth <- function(z) {
    stopifnot(is.numeric(z) || is.complex(z))
    atanh(1/z)
}

# area cosecans
acsch <- function(z) {
    stopifnot(is.numeric(z) || is.complex(z))
    asinh(1/z)
}

# area secans
asech <- function(z) {
    stopifnot(is.numeric(z) || is.complex(z))
    acosh(1/z)
}
