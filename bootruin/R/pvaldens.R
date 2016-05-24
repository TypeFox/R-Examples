pvaldens <- function(x, bw, rho, method = c("jh", "chen")) {
    stopifnot(is.numeric(x))

    x <- as.vector(x)
    x <- x[is.finite(x)]

    method <- match.arg(method)
    switch(method,
        jh = {
            if (missing(rho)) {
                rho <- ifelse(missing(bw), 0.9, 1.0 - bw^2.0)
            } else {
                if (any(rho <= 0.0, rho >= 1.0)) {
                    stop("'rho' must have a value between 0 and 1.")
                }
            }

            K.jh <- function(u, v){
                if(any(u <= 0.0, v <= 0.0, u >= 1.0, v >= 1.0)){
                    return(0.0)
                } else {
                    1.0 / sqrt(1.0 - rho^2.0) *
                    exp(-rho / (2.0 * (1.0 - rho^2.0)) *
                    (rho * qnorm(u)^2.0 + rho * qnorm(v)^2.0 -
                    2.0 * qnorm(u) * qnorm(v)))
                }
            }

            return(
                function(y) {
                    rowMeans(outer(y, x, Vectorize(K.jh)))
                }
            )
        },
        chen = {
            if (missing(bw)) {
                bw <- 0.01
            }

            if (!missing(rho)) {
                warning(paste("The parameter 'rho' is not supported by the method '", method, "'.", sep = ""))
            }

            return(
                function(y) {
                    colMeans(outer(x, y, function(x, y) {
                                   dbeta(y, x / bw + 1.0, (1.0 - x) / bw + 1.0)
                             }),
                             na.rm = TRUE)
                }
            )
        }
    )
}
