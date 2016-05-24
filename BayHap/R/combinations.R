'combinations'<-function (n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE) 
{
    if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n%%1) != 
        0) 
        stop("bad value of n")
    if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r%%1) != 
        0) 
        stop("bad value of r")
    if (!is.atomic(v) || length(v) < n) 
        stop("v is either non-atomic or too short")
    if ((r > n) & repeats.allowed == FALSE) 
        stop("r > n and repeats.allowed=FALSE")
    if (set) {
        v <- unique(sort(v))
        if (length(v) < n) 
            stop("too few different elements")
    }
    v0 <- vector(mode(v), 0)
    if (repeats.allowed) 
        sub <- function(n, r, v) {
            if (r == 0) 
                v0
            else if (r == 1) 
                matrix(v, n, 1)
            else if (n == 1) 
                matrix(v, 1, r)
            else rbind(cbind(v[1], Recall(n, r - 1, v)), Recall(n - 
                1, r, v[-1]))
        }
    else sub <- function(n, r, v) {
        if (r == 0) 
            v0
        else if (r == 1) 
            matrix(v, n, 1)
        else if (r == n) 
            matrix(v, 1, n)
        else rbind(cbind(v[1], Recall(n - 1, r - 1, v[-1])), 
            Recall(n - 1, r, v[-1]))
    }
    sub(n, r, v[1:n])
}

