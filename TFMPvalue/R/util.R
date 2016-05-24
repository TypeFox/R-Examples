### Typical 'prior.params' vector: c(A=0.25, C=0.25, G=0.25, T=0.25)
### This is taken from Biostrings package matchPWM.R.
### Just to get rid of the node during the build.
DNA_BASES <- c("A", "C", "G", "T")

normargPriorParams <- function(prior.params)
{
    if (!is.numeric(prior.params))
        stop("'prior.params' must be a numeric vector")
    if (length(prior.params) != length(DNA_BASES) ||
        !setequal(names(prior.params), DNA_BASES))
        stop("'prior.params' elements must be named A, C, G and T")
    ## Re-order the elements.
    prior.params <- prior.params[DNA_BASES]
    if (any(is.na(prior.params)) || any(prior.params < 0))
        stop("'prior.params' contains NAs and/or negative values")
    prior.params
}

normargMat = function(x){
    if (is.null(rownames(x)))
        stop("invalid Matrix 'mat': no row names")
    if (!all(DNA_BASES %in% rownames(x)))
        stop("invalid Matrix 'mat': row names must contain A, C, G and T")
    if (any(duplicated(rownames(x))))
        stop("invalid Matrix 'mat': duplicated row names")
    if (ncol(x) == 0L)
        stop("invalid Matrix 'mat': no columns")
    if (any(is.na(x)) || any(x < 0L))
        stop("invalid Matrix 'mat': values cannot be NA or negative")
    if (any(x[!(rownames(x) %in% DNA_BASES), ] != 0L))
        stop("invalid Matrix 'mat': IUPAC ambiguity letters are represented")
    x <- x[DNA_BASES, , drop = FALSE]
    x
}

