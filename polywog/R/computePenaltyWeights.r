##
## Compute the penalty weights for a polywog model (NULL if fitting via SCAD)
##
computePenaltyWeights <- function(X,
                                  y,
                                  weights,
                                  polyTerms,
                                  lmcoef,
                                  method,
                                  penwt.method,
                                  family,
                                  unpenalized)
{
    if (method == "scad") {
        ## No weights when using SCAD
        ans <- NULL
    } else if (family == "gaussian" || penwt.method == "lm") {
        ## Inverse linear model coefficients for continuous outcomes or when
        ## specified for a binary model
        ans <- 1 / abs(lmcoef)[-1]
        ans[unpenalized] <- 0
    } else {
        ## Inverse logit coefficients for a binary model
        ans <- penaltyWeightsBinary(X, y, polyTerms, weights)
        ans[unpenalized] <- 0
    }

    ans
}
