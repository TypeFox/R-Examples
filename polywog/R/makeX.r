##
## Assemble the 1-dimensional (i.e., non-polynomial-expanded) model matrix
##
## ARGUMENTS:
##   formula: model formula (of class "Formula")
##   mf: model frame
##
## RETURN:
##   model matrix (no intercept), with attributes "k_expand" telling how many
##   columns should be included in the polynomial expansion and "k_lin"
##   telling how many are included linearly
##
makeX <- function(formula, mf)
{
    ## Ensure no intercepts included
    formula <- removeIntercepts(formula, mf)

    ## Extract the model matrix from the model frame
    X <- model.matrix(formula, data = mf, rhs = 1)
    k_expand <- ncol(X)
    if (length(formula)[2] > 1)
        X <- cbind(X, model.matrix(formula, data = mf, rhs = 2))
    k_lin <- ncol(X) - k_expand
    binary_cols <- apply(X, 2, function(a) all(a %in% 0:1))

    structure(X,
              k_expand = k_expand,
              k_lin = k_lin,
              binary_cols = binary_cols)
}
