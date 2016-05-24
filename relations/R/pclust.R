relation_pclust <-
function(x, k, method, m = 1, weights = 1, control = list())
{
    x <- as.relation_ensemble(x)

    ## Handle 'method'.
    if(is.character(method))
        method <- get_relation_consensus_method(method)
    if(inherits(method, "relation_consensus_method")) {
        ## Not all consensus methods (e.g., Borda) give central
        ## relations.
        if(is.null(md <- method$dissimilarity)
           || is.null(me <- method$exponent))
            stop("Consensus method does not compute central relations.")
        family <-
            clue::pclust_family(D = function(x, y = NULL)
                                relation_dissimilarity(x, y, md),
                                C = method$definition,
                                e = me)
    }
    else if(inherits(method, "pclust_family")) {
        ## User-defined, maybe.
        ## <FIXME>
        ## Eventually, this should be deprecated in favor of registering
        ## consensus methods giving central relations.
        family <- method
        ## </FIXME>
    }
    else
        stop("Invalid 'method' argument.")

    out <- clue::pclust(x, k, family, m, weights, control)

    ## Massage the results a bit.
    dissimilarities <- as.matrix(family$D(x) ^ family$e)
    ## <FIXME>
    ## Component 'prototypes' should really be a relation ensemble.
    ## But don't we get this automatically?
    ##   out$prototypes <- relation_ensemble(list = out$prototypes)
    ## </FIXME>
    out$call <- match.call()
    out <- .structure(c(out,
                        list(silhouette =
                             silhouette(out$cluster,
                                        dmatrix = dissimilarities),
                             validity =
                             clue::cl_validity(out$membership,
                                               dissimilarities))),
                      class =
                      unique(c("relation_pclust", class(out))))

    clue::as.cl_partition(out)
}

print.relation_pclust <-
function(x, ...)
{
    txt <- if(x$m == 1)
        gettextf("A hard partition of a relation ensemble with %d elements into %d classes.",
                 clue::n_of_objects(x), clue::n_of_classes(x))
    else
        gettextf("A soft partition (degree m = %f) of a relation ensemble with %d elements into %d classes.",
                 x$m, clue::n_of_objects(x), clue::n_of_classes(x))
    writeLines(strwrap(txt))
    NextMethod("print", x, header = FALSE)
    print(x$validity, ...)
    invisible(x)
}
