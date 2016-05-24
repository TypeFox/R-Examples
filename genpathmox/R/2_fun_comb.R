#' @title Enumerating the possible combinations of a specified size from the elements of a vector
#' @details
#' Internal function. \code{comb} is called by \code{partition}.
#' @param n Size of the source vector
#' @param r Size of the target vectors
#' @param v Source vector. Defaults to 1:n
#' @param set Logical flag indicating whether duplicates should be removed from the source vector v
#' @param repeats.allowed Logical flag indicating whether the constructed vectors may include duplicated values.
#' @return  matrix
#' @keywords internal
#' @export

comb	<-	function (n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE)
{
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
