deffK <- function(w){
    if (any(w <= 0))
        warning("Some weights are less than or equal to 0.\n")
    n <- length(w)
    1 + sum((w - mean(w))^2)/n/mean(w)^2
}
