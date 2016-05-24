#' Disattenuate a correlation matrix using an estimate of the component reliabilities
#' 
#' @param r_mat A correlation matrix
#' @param rel_vec A vector or reliabilities.
#' @return A reliabated (disattenauted) correlation matrix.
#' @author Allen Goebl and Jeff Jones
#' @examples
#' r_mat <- matrix(c(1.00, 0.25, 0.30, 
#'                   0.25, 1.00, 0.50, 
#'                   0.30, 0.50, 1.00), 3, 3)
#' rel   <- c(.70, .64, .81)
#' reliabate(r_mat = r_mat, rel_vec = rel)
#' @export
reliabate <- function (r_mat, rel_vec) {
    out <- diag(1/sqrt(rel_vec)) %*% r_mat %*% diag(1/sqrt(rel_vec))
    diag(out) <- 1
    return(out)
}
