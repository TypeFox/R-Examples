"QUnif" <-
function (n, min = 0, max = 1, n.min = 1, p, leap = 1) 
{
    stopifnot(1 <= (p <- as.integer(p)), 1 <= (n <- as.integer(n)), 
        1 <= (leap <- as.integer(leap)), 1 <= (n.min <- as.integer(n.min)))
    stopifnot((n.max <- n.min + (n - 1:1) * leap) < .Machine$integer.max)
    pr. <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 
        43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 
        103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 
        163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 
        227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 
        281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 
        353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 
        421, 431, 433, 439, 443, 449, 457)
    if (length(pr.) < p) 
        stop("primes not yet available for p=", p)
    pr <- pr.[1:p]
    if (leap > 1 && any(leap == pr) && length(pr.) >= p + 1) 
        pr <- c(pr[leap != pr], pr.[p + 1])
    stopifnot(length(max) == p || length(max) == 1, length(min) == 
        p || length(min) == 1)
    max <- rep.int(max, p)
    min <- rep.int(min, p)
    dU <- max - min
    r <- matrix(0, n, p)
    for (j in 1:p) r[, j] <- min[j] + dU[j] * sHalton(n.max, 
        n.min, base = pr[j], leap = leap)
    r
}

