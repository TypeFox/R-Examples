base <-
function (n, base = 10, num.digits = max(0, floor(log(n, base))) + 
    1) 
{
    if (!is.vector(n, mode = "numeric") || is.factor(n) || length(n) != 
        1 || n != trunc(n) || n < 0) 
        stop("'n' must be a non-negative integer")
    if (!is.vector(base, mode = "numeric") || is.factor(base) || 
        length(base) != 1 || base != trunc(base) || base <= 1) 
        stop("'base' must be positive integer greater than 1")
    if (n == 0) {
        if (!is.vector(num.digits, mode = "numeric") || is.factor(num.digits) || 
            length(num.digits) != 1 || num.digits != trunc(num.digits) || 
            num.digits < 1) 
            stop("'num.digits' must be a positive integer")
        vec <- numeric(num.digits)
    }
    else {
        min.num.digits <- floor(log(n, base)) + 1
        if (!is.vector(num.digits, mode = "numeric") || is.factor(num.digits) || 
            length(num.digits) != 1 || num.digits != trunc(num.digits) || 
            num.digits < min.num.digits) 
            stop(paste("'num.digits' must be a positive integer", 
                "greater than or equal to", min.num.digits))
        vec <- numeric(num.digits)
        vec[(num.digits - min.num.digits) + (1:min.num.digits)] <- (n%/%(base^((min.num.digits - 
            1):0)))%%base
    }
    vec
}
