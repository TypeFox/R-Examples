##
##  n o r m . R  Vector Norm
##


Norm <- function(x, p=2) {
    stopifnot(is.numeric(x) || is.complex(x),
              is.numeric(p), length(p) == 1)

    if (p > -Inf && p < Inf) sum(abs(x)^p)^(1/p)
    else if (p ==  Inf) max(abs(x))
    else if (p == -Inf) min(abs(x))
    else return(NULL)
}
