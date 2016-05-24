##
## Compute the polynomial expansion of a matrix, according to the specified
## `polyTerms` object
##
expandMatrix <- function(X, polyTerms, intercept = FALSE)
{
    ans <- computeExpandMatrix(X = X, poly_terms = polyTerms)
    rownames(ans) <- rownames(X)
    colnames(ans) <- rownames(polyTerms)
    if (intercept)
        ans <- cbind("(Intercept)" = 1L, ans)
    ans
}
