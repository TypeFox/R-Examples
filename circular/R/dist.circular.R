#############################################################
#                                                           #
#   dist.circular function                                  #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: October, 12, 2009                                 #
#   Version: 0.1-2                                          #
#                                                           #
#   Copyright (C) 2007 Claudio Agostinelli                  #
#                                                           #
#############################################################

dist.circular <- function (x, method = "correlation", diag = FALSE, upper = FALSE) {
    x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
    attr(x, "class") <- attr(x, "circularp") <-  NULL
    if (!is.na(pmatch(method, "correlation"))) 
        method <- "correlation"
    METHODS <- c("correlation", "angularseparation", "chord", "geodesic")
    method <- pmatch(method, METHODS)
    if (is.na(method)) 
        stop("invalid distance method")
    if (method == -1) 
        stop("ambiguous distance method")
    N <- nrow(x <- as.matrix(x))
    d <- .C("R_distance", x = as.double(x), nr = N, nc = ncol(x), 
        d = double(N * (N - 1)/2), diag = as.integer(FALSE), 
        method = as.integer(method), DUP = FALSE, 
        NAOK = TRUE, PACKAGE = "circular")$d
    attr(d, "Size") <- N
    attr(d, "Labels") <- dimnames(x)[[1]]
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- METHODS[method]
#    if (method == 6) 
#        attr(d, "p") <- p
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    return(d)
}
