##
## Compute the "polyTerms" matrix
##
## ARGUMENTS:
##   degree: degree of polynomial expansion
##   k_expand: number of terms to include in the expansion
##   k_lin: number of terms to include linearly
##   binary_cols: logical vector indicating which columns of the raw input
##     matrix contain binary variables
##   names.: optional character vector containing the name of each raw input
##     variable
##
## RETURN:
##   Matrix where each column is a raw model term, each row is a term of the
##   polynomial expansion, and each entry ij is the degree of raw term i in
##   expansion term j
##
makePolyTerms <- function(degree, k_expand, k_lin, binary_cols, names. = NULL)
{
    ## Call C++ routine to compute the full set of potential polynomial terms
    ans <- computePolyTerms(degree, k_expand, k_lin)

    ## Eliminate duplicates
    ##
    ## This is done before combining the elements of ans, since a degree-k row
    ## can't be a duplicate of a degree-j row if k != j
    ans <- lapply(ans, function(x) x[!duplicated(x), , drop = FALSE])

    ## Combine into a single matrix
    ans <- do.call(rbind, ans)

    ## Remove any expanded term that contains a square power (or greater) of a
    ## binary variable, since it will be collinear with a lower-order term
    if (any(binary_cols)) {
        any_higher_binary <- rowSums(ans[, binary_cols, drop = FALSE] > 1)
        ans <- ans[!any_higher_binary, , drop = FALSE]
    }

    ## Add row and column names (if requested)
    if (!is.null(names.)) {
        ## Column names
        colnames(ans) <- names.

        ## Row names
        termNames <- apply(ans, 1, function(x) {
            deg <- as.character(x)
            deg <- ifelse(x == 1,
                          "",
                          paste("^", deg, sep = ""))
            deg <- paste(names., deg, sep = "")
            deg <- paste(deg[x > 0], collapse = ".")
            deg
        })
        rownames(ans) <- termNames
    }

    ans
}
