lspline <- function(x, df = NULL, knots = NULL) {
        x <- as.numeric(x)
        
        if(is.null(knots) && is.null(df))
                return(x)
        nax <- is.na(x)
        
        if(hasNA <- any(nax))
                x <- x[!nax]
        if(is.null(knots))
                knots <- quantile(x, seq(0, 1, len = df + 2 - 1))
        B <- splineDesign(knots, x, ord = 1)

        if(!isTRUE(NROW(B) == length(x)))
                stop("problem with spline design matrix")
        basis <- B * x
        
        if(hasNA) {
                mat <- matrix(NA, nrow = length(nax), ncol = ncol(basis))
                mat[!nax, ] <- basis
                basis <- mat
        }
        structure(basis, knots = knots, class = "basis")
}
