#' Create a response style
#' 
#' Creates a response style by cutting up a quadratic monotone spline.
#' 
#' @param alpha vector of spline coefficients
#' @param nr.scale number of rating categories; numeric
#' @param tvec knots for spline functions
#' @param xvec evaluation points for basis functions
#' @param scale logical; scale or not
#' @author Pieter C. Schoonees
#' @keywords splines
#' @export create.rs
create.rs <-
function (alpha = matrix(c(1, 2, 1), nrow = 1), nr.scale = 7, 
          tvec = c(0, 0.5, 1), xvec = 0:nr.scale/nr.scale, scale = TRUE) 
{
    if(!is.matrix(alpha)) alpha <- as.matrix(alpha)
    if(ncol(alpha) == 1) alpha <- t(alpha)
    if(ncol(alpha) == 4) {
      cat("Dropping first column of alpha: assuming an intercept...\n\n")
      alpha <- alpha[,-1]
    }
    if(ncol(alpha) >= 5) stop("Incorrect number of columns\n")
    create.rs.vec <- function (x){
      M <- ispline(xvec = xvec, tvec = tvec, 
                   intercept = FALSE)
      out <- as.numeric(M %*% matrix(x, ncol = 1)/
                          ifelse(scale, sum(x), 1))
      out
    }
    out <- t(apply(alpha, 1, function(x) create.rs.vec(x)))
    if(nrow(alpha) == 1) out <- as.numeric(out)
    rownames(out) <- rownames(alpha)
    out
}
