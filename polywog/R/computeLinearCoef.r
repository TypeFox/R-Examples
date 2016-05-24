##
## Compute linear model coefficients and pivot columns via the QR
## decomposition of the polynomial-expanded model matrix
##
computeLinearCoef <- function(X,
                              y,
                              polyTerms,
                              weights,
                              allowRankDeficient = TRUE)
{
    X.expand <- expandMatrix(X = X,
                             polyTerms = polyTerms,
                             intercept = TRUE)
    qx <- qr(sqrt(weights) * X.expand)

    ## Check for rank deficiency and either throw an error or find the pivot
    ## columns
    if (!allowRankDeficient && qx$rank < ncol(X.expand)) {
        stop("Model matrix is rank-deficient")
    } else {
        pivot <- qx$pivot[seq_len(qx$rank)]
        lmcoef <- qr.coef(qx, sqrt(weights) * y)[pivot]
        pivot <- pivot[-1] - 1L  # Remove intercept from pivot columns
    }

    structure(lmcoef, pivot = pivot)
}
