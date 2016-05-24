### * relation

relation <-
function(domain = NULL, incidence = NULL, graph = NULL, charfun = NULL)
{
    if(sum(is.null(incidence), is.null(graph), is.null(charfun)) != 2L)
        stop("Need exactly one of 'incidence', 'graph', and 'charfun'.")

    if(!is.null(domain)) {
        ## Be nice first ...
        if(!is.list(domain) || is.cset(domain))
            domain <- list(X = domain)
        ## ... and then check.
        if(!.is_valid_relation_domain(domain))
            stop("Invalid relation domain.")
    }

    if(!is.null(incidence)) {
        incidence <- as.array(incidence)
        if(!.is_valid_relation_incidence(incidence))
            stop("Invalid relation incidence.")
        size <- dim(incidence)
        if(!is.null(domain)) {
            ## Be nice first ...
            domain <- rep(domain, length.out = length(size))
            ## ... and then check.
            if(any(size != sapply(domain, length)))
                stop("Relation size mismatch between domain and incidence.")
        }
        return(.make_relation_from_domain_and_incidence(domain, incidence))
    }

    if(!is.null(graph)) {
        if (is.gset(graph) &&
            (gset_is_multiset(graph) || gset_is_fuzzy_multiset(graph)))
            stop("Only crisp or fuzzy sets allowed.")
        G <- .make_relation_graph_components(graph)
        ## Be nice and recycle domain (useful for endorelations).
        if (!is.null(domain) && (length(G) > 0L))
            domain <- rep(domain, length.out = length(G))
        return(.make_relation_from_domain_and_graph_components(domain, G))
    }

    if(!is.null(charfun)) {
        if(is.null(domain))
            stop("Need domain along with characteristic function.")
        ## No recycling here, as we really do not know the arity of a
        ## function (nor is this a well-defined notion).
        I <- array(do.call(mapply,
                           c(list(charfun),
                             .cartesian_product(lapply(domain, as.list)))),
                   dim = sapply(domain, length))
        return(.make_relation_from_domain_and_incidence(domain, I))
    }
}

homorelation <-
function(domain = NULL, incidence = NULL, graph = NULL, charfun = NULL)
{
    if(sum(is.null(incidence), is.null(graph), is.null(charfun)) != 2L)
        stop("Need exactly one of 'incidence', 'graph', and 'charfun'.")

    if(!is.null(domain)) {
        arity <- length(domain)

        ## merge domain labels
        domain <- do.call(cset_union, domain)

        ## recycle domain for charfun-generators
        if (!is.null(charfun))
            domain <- rep.int(list(domain), arity)
    } else {
        if(!is.null(graph)) {
            ## merge domain labels
            domain <- sort(unique(unlist(graph)))

            ## for data frame graphs, fix domain names
            if(!is.null(graph) && is.data.frame(graph))
                names(graph) <- rep.int("X", ncol(graph))
        } else if(!is.null(incidence)) {
            dn <- dimnames(incidence)

            ## merge domain labels taken from array dimnames
            domain <- unique(unlist(dn))

            ## match merged domain against actual labels
            ind <- Map(match, list(domain), dn)

            ## span target array using indices
            incidence <- do.call("[", c(list(incidence), ind))

            ## replace all NAs produced with 0
            incidence[is.na(incidence)] <- 0
        }
    }

    R <- relation(domain = domain, incidence = incidence,
                  graph = graph, charfun = charfun)
    .set_property(R, "is_homogeneous", TRUE)
}

endorelation <-
function(domain = NULL, incidence = NULL, graph = NULL, charfun = NULL)
{
    if ((!is.null(incidence) && length(dim(incidence)) != 2L)
        || (!is.null(graph)
            && is.data.frame(graph) && ncol(graph) != 2L)
        || (!is.null(graph)
            && is.set(graph)
            && any(sapply(graph, length) != 2L))
        || (!is.null(charfun) && length(domain) != 2L ))
        stop("Relation is not binary.")
    R <- homorelation(domain = domain, incidence = incidence,
                      graph = graph, charfun = charfun)
    .set_property(R, "is_endorelation", TRUE)
}

### * is.relation

is.relation <-
function(x)
    inherits(x, "relation")

### * as.relation

as.relation <-
function(x, ...)
    UseMethod("as.relation")

as.relation.default <-
function(x, ...)
    stop("Method not implemented.")

## Obviously.
as.relation.relation <-
function(x, ...) x

## Logical vectors are taken as unary relations (predicates).
as.relation.logical <-
function(x, ...)
{
    D <- if(!is.null(nms <- names(x)) && !any(duplicated(nms)))
        list(nms)
    else
        NULL
    I <- as.array(as.integer(x))
    meta <- list(is_endorelation = FALSE,
                 is_complete = all(is.finite(x)))
    .make_relation_from_domain_and_incidence(D, I, meta)
}

## Numeric vectors and ordered factors are taken as order relations.
as.relation.numeric <-
function(x, ...)
{
    D <- if(!is.null(nms <- names(x)) && !any(duplicated(nms)))
        list(nms, nms)
    else if(!any(duplicated(x)))
        rep.int(list(x), 2L)
    else
        NULL
    I <- outer(x, x, "<=")
    meta <- if(any(is.na(x)))
        list(is_endorelation = TRUE,
             is_complete = NA,
             is_reflexive = NA,
             is_antisymmetric = NA,
             is_transitive = NA)
    else
        list(is_endorelation = TRUE,
             is_complete = TRUE,
             is_reflexive = TRUE,
             is_antisymmetric = !any(duplicated(x)),
             is_transitive = TRUE)
    .make_relation_from_domain_and_incidence(D, I, meta)
}
as.relation.integer <- as.relation.numeric

## This is almost identical to as.relation.numeric, but when generating
## domains from unique values we use these as is in the numeric case,
## but use as.character() of these for (ordered) factors.
as.relation.ordered <-
function(x, ...)
{
    D <- if(!is.null(nms <- names(x)) && !any(duplicated(nms)))
        list(nms, nms)
    else if(!any(duplicated(x)))
        rep.int(list(as.character(x)), 2L)
    else
        NULL
    I <- outer(x, x, "<=")
    meta <- if(any(is.na(x)))
        list(is_endorelation = TRUE,
             is_complete = NA,
             is_reflexive = NA,
             is_antisymmetric = NA,
             is_transitive = NA)
    else
        list(is_endorelation = TRUE,
             is_complete = TRUE,
             is_reflexive = TRUE,
             is_antisymmetric = !any(duplicated(x)),
             is_transitive = TRUE)
    .make_relation_from_domain_and_incidence(D, I, meta)
}

## Unordered factors are taken as equivalence relations.
as.relation.factor <-
function(x, ...)
{
    D <- if(!is.null(nms <- names(x)) && !any(duplicated(nms)))
        list(nms, nms)
    else if(!any(duplicated(x)))
        rep.int(list(as.character(x)), 2L)
    else
        NULL
    I <- outer(x, x, "==")
    meta <- if(any(is.na(x)))
        list(is_endorelation = TRUE,
             is_reflexive = NA,
             is_symmetric = NA,
             is_transitive = NA)
    else
        list(is_endorelation = TRUE,
             is_complete = isTRUE(all(x == x[1L])),
             is_reflexive = TRUE,
             is_symmetric = TRUE,
             is_transitive = TRUE)
    .make_relation_from_domain_and_incidence(D, I, meta)
}

as.relation.character <-
function(x, ...)
    as.relation(factor(x))

as.relation.set <-
function(x, ...)
    relation(graph = x)

as.relation.gset <-
function(x, ...)
    relation(graph = x)

## Matrices and arrays are taken as incidences of relations, provided
## that they are feasible.
as.relation.matrix <-
function(x, ...)
{
    if(!.is_valid_relation_incidence(x))
        stop("Invalid relation incidence.")
    meta <- list(is_endorelation =
                 .relation_is_endorelation_using_incidence(x))
    .make_relation_from_domain_and_incidence(dimnames(x), x, meta)
}
as.relation.array <- as.relation.matrix

## Data frames are taken as relation graph components.
as.relation.data.frame <-
function(x, ...)
{
    ## Get the domain.
    D <- lapply(x, unique)
    names(D) <- names(x)

    ## Get the incidences.
    I <- .make_incidence_from_domain_and_graph_components(D, x)

    ## use membership information, if any
    M <- attr(x, "memberships")
    if (!is.null(M))
        I[I > 0] <- M

    ## And put things together.
    .make_relation_from_domain_and_incidence(D, I)
}

## Package clue: cl_partition objects.
## <FIXME>
## Of course, CLUE has a notion of soft ("fuzzy") partitions in the
## Ruspini sense, but it is not clear how these correspond (should be
## mapped) to fuzzy equivalence relations.
## </FIXME>
as.relation.cl_partition <-
function(x, ...)
    as.relation(factor(clue::cl_class_ids(x)))

## Package seriation: ser_permutation objects.
as.relation.ser_permutation <-
function(x, ...)
{
    o <- seriation::get_order(x)
    oo <- order(o)
    D <- if(!is.null(nms <- names(o)[oo]) && !any(duplicated(nms)))
        list(nms, nms)
    else
        NULL
    I <- outer(oo, oo, "<=")
    meta <- .relation_meta_db[["L"]]
    .make_relation_from_domain_and_incidence(D, I, meta)
}


### * Relation methods

### ** [.relation

`[.relation` <-
function(x, ...)
{
    l <- match.call()[- (1 : 2)]
    D <- .domain(x)
    I <- relation_incidence(x)
    d <- dim(I)
    dn <- dimnames(I)

    if (length(l) != length(D))
        stop("Wrong number of arguments.")

    ## complete empty indices
    l <- lapply(seq_along(l), function(i) {
        ## check missings
        if (identical(l[[i]], alist(, )[[1L]]))
            return(seq_len(d[i]))
        ## retrieve parameter
        e <- eval(l[[i]])
        ## integer arg -> use as index
        if (is.numeric(e) && (e == as.integer(e)))
            e
        ## character arg -> match against dimnames labels
        else if (is.character(e))
            match(e, dn[[i]])
        ## else: match against domain elements
        else
            .exact_match(list(e), D[[i]])
    })

    ## subset domain
    D <- Map(.set_subset, D, l)

    ## subset incidence matrix using standard method
    I <- do.call("[", c(list(I), l, drop = FALSE))

    ## return new relation
    .make_relation_from_domain_and_incidence(D, I)
}

### ** all.equal.relation

all.equal.relation <-
function(target, current, check.attributes = TRUE, ...)
{
    ## Note that we really do not know what 'current' is.  So we compare
    ## classes before anything else.
    if(data.class(target) != data.class(current)) {
        ## Common msg style, but i18ned.
        return(gettextf("target is %s, current is %s",
                        data.class(target), data.class(current)))
    }

    msg <- if(check.attributes)
        attr.all.equal(relation_properties(target),
                       relation_properties(current), ...)

    ## Compare arities, then sizes, then domains.
    D_t <- relation_domain(target)
    D_c <- relation_domain(current)
    a_t <- length(D_t)
    a_c <- length(D_c)
    if(a_t != a_c)
        return(c(msg,
                 gettextf("Relation arities (%d, %d) differ.",
                          a_t, a_c)))
    s_t <- sapply(D_t, length)
    s_c <- sapply(D_c, length)
    if(!identical(s_t, s_c))
        return(c(msg,
                 gettextf("Relation sizes (%s, %s) differ.",
                          paste(s_t, collapse = "/"),
                          paste(s_c, collapse = "/"),
                          domain = NA)))
    if(!.domain_is_equal(D_t, D_c)) {
        ## Maybe use all.equal.set eventually.
        return(c(msg,
                 gettextf("Relation domains differ in elements.")))
    }
    aei <- all.equal(relation_incidence(target),
                     relation_incidence(current))
    if(!identical(aei, TRUE))
        aei <- c("Relation incidences differ:", aei)
    c(msg, aei)
}

### ** as.data.frame.relation

as.data.frame.relation <-
function(x, row.names = NULL, ...)
{
    ## Get the "raw" graph components.
    out <- .make_relation_graph_components(x)
    M <- attr(out, "memberships")
    names(out) <-
        .make_domain_names_from_relation_graph_components(out,
                                                          relation_is_endorelation(x))
    ## Flatten.
    out <- lapply(out, unlist, recursive = FALSE)

    ## And put into "some kind of" data frame.
    .structure(.make_data_frame_from_list(out, row.names),
               memberships = M)
}

### ** as.tuple.relation

as.tuple.relation <-
function(x)
{
    D <- as.tuple(relation_domain(x))
    G <- as.set(relation_graph(x))
    if(is.null(names(D)))
        names(D) <- names(G)
    names(G) <- NULL
    pair(Domain = D, Graph = G)
}

### * cut.relation

cut.relation <-
function(x, level = 1, ...)
    .make_relation_from_domain_and_incidence(.domain(x),
                                             .incidence(x) >= level)

### * dim.relation

dim.relation <- relation_size

### ** print.relation

print.relation <-
function(x, ...)
{
    a <- .arity(x)
    s <- paste(.size(x), collapse = " x ")
    if(identical(relation_is_crisp(x), FALSE)) {
        if(a == 1L)
            writeLines(gettextf("A unary fuzzy relation of size %s.", s))
        else if(a == 2L)
            writeLines(gettextf("A binary fuzzy relation of size %s.", s))
        else
            writeLines(gettextf("A %d-ary fuzzy relation of size %s.", a, s))
    } else {
        if(a == 1L)
            writeLines(gettextf("A unary relation of size %s.", s))
        else if(a == 2L)
            writeLines(gettextf("A binary relation of size %s.", s))
        else
            writeLines(gettextf("A %d-ary relation of size %s.", a, s))
    }
    invisible(x)
}

summary.relation <-
function(object, ...)
{
    .structure(.check_all_predicates(object, ...),
               class = "summary.relation")
}

print.summary.relation <-
function(x, ...)
{
    print(unclass(x))
}

### * Group methods and related.

## Here is what we do.

## * Comparisons are obvious.
## * We use & and | for intersection and union.
## * We use min/max for meet and join (which of course is the same as
##   the above).
## * We use * for the composition and unary ! for the converse.
## * Finally, t() is used for the inverse.

Summary.relation <-
function(..., na.rm = FALSE)
{
    ok <- switch(.Generic, max = , min = , range = TRUE, FALSE)
    if(!ok)
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class),
             domain = NA)
    args <- list(...)
    x <- relation_ensemble(list = args)
    switch(.Generic,
           "min" = .relation_meet(x),
           "max" = .relation_join(x),
           "range" = {
               relation_ensemble(min = .relation_meet(x),
                                 max = .relation_join(x))
           })
}


Ops.relation <-
function(e1, e2)
{
    if(nargs() == 1L) {
        if(!(as.character(.Generic) %in% "!"))
            stop(gettextf("Unary '%s' not defined for \"%s\" objects.",
                          .Generic, .Class),
                 domain = NA)
        return(relation_complement(e1))
    }

    ## In addition to comparisons, we support & | * + - %/% %% and ^.
    if(!(as.character(.Generic)
         %in% c("<", "<=", ">", ">=", "==", "!=",
                "&", "|", "*", "+", "-", "%/%", "%%", "^")))
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class),
             domain = NA)

    switch(.Generic,
           "+" = return(relation_union(e1, e2)),
           "-" = return(relation_complement(e1, e2)),
           "%/%" = return(relation_division(e1, e2)),
           "%%" = return(relation_remainder(e1, e2))
           )

    D1 <- .domain(e1)
    I1 <- .incidence(e1)

    if(as.character(.Generic == "^")) {
        ## Allow for nonnegative integer powers of endorelations.
        if(!relation_is_endorelation(e1))
            stop("Power only defined for endorelations.")
        if((length(e2) != 1L) || (e2 < 0) || (e2 != round(e2)))
            stop("Power only defined for nonnegative integer exponents.")
        if(e2 == 0) {
            ## Return the equality relation: maybe encapsulate this, or
            ## at least add metadata?
            I <- diag(nrow = .size(e1)[[1L]])
            return(.make_relation_from_domain_and_incidence(D1, I))
        } else {
            I <- Reduce("%*%", rep.int(I1, e2))
            return(.make_relation_from_domain_and_incidence(D1, I))
        }
    }

    D2 <- .domain(e2)
    I2 <- .incidence(e2)

    ## Composition (*) is only defined for binary relations with
    ## matching 2nd/1st domain elements.
    if(as.character(.Generic) == "*") {
        if((length(D1) != 2L)
           || (length(D2) != 2L)
           || !set_is_equal(D1[[2L]], D2[[1L]]))
            stop("Composition of given relations not defined.")
        ## When composing indicidences, need the same *internal* order
        ## for D2[[1L]] as for D1[[2L]].
        if(isTRUE(relation_is_crisp(e1)) && isTRUE(relation_is_crisp(e2)))
            I <- ((I1 %*% I2) > 0)
        else {
            n <- ncol(I1)               # Same as nrow(I2).
            I <- matrix(0, nrow = nrow(I1), ncol = ncol(I2))
            for(j in seq_len(n))
                I <- pmax(I, outer(I1[, j], I2[j, ], .T.))
        }
        ## The composition has domain (D1[[1L]], D2[[2L]]) (and
        ## appropriate names).  Information about auto-generation of
        ## domains is currently ignored.
        D <- list(D1[[1L]], D2[[2L]])
        if(!is.null(nms <- names(D1)))
            names(D)[1L] <- nms[1L]
        if(!is.null(nms <- names(D2)))
            names(D)[2L] <- nms[2L]
        return(.make_relation_from_domain_and_incidence(D, I))
    }

    ## In the remaining cases, the relations must have equal domains.
    if(!.domain_is_equal(D1, D2))
        stop("Relations need equal domains.")
    switch(.Generic,
           "<=" = all(I1 <= I2),
           "<"  = all(I1 <= I2) && any(I1 < I2),
           ">=" = all(I1 >= I2),
           ">"  = all(I1 >= I2) && any(I1 > I2),
           "==" = all(I1 == I2),
           "!=" = any(I1 != I2),
           "&" = .make_relation_from_domain_and_incidence(D1, .T.(I1, I2)),
           "|" = .make_relation_from_domain_and_incidence(D1, .S.(I1, I2))
           )
}

rev.relation <-
t.relation <-
function(x)
{
    if(!relation_is_binary(x))
        stop("Can only compute inverses of binary relations.")
    .make_relation_from_domain_and_incidence(rev(.domain(x)),
                                             t(.incidence(x)))
}

### * Relation representations

### ** .make_relation_by_domain_and_incidence

.make_relation_by_domain_and_incidence <-
function(D, I)
{
    ## Generate a valid domain if needed.
    if(!.is_valid_relation_domain(D)) {
        D <- dimnames(I)
        if(!.is_valid_relation_domain(D)) {
            D <- lapply(dim(I), function(s) as.character(seq_len(s)))
            ## Just to make a point ...
            attr(D, "auto") <- TRUE
        }
    }

    size <- dim(I)
    ## Strip incidence of all attributes but dim.
    I <- .structure(c(I), dim = size)

    ## Now canonicalize by turning all domain components into sets, and
    ## reordering the incidences accordingly.  Note that components
    ## which are already sets are already in the canonical order.
    ind <- sapply(D, is.cset)
    if(!all(ind)) {
        ## If all components are sets there is nothing we need to do.
        pos <- vector("list", length = length(D))
        if(any(ind)) pos[ind] <- lapply(D[ind], seq_along)
        ind <- !ind
        sets_with_order <- lapply(D[ind], make_set_with_order)
        ## Turn all non-set domain components into sets.
        D[ind] <- lapply(sets_with_order, `[[`, "set")
        ## Reorder incidences accordingly.
        pos[ind] <- lapply(sets_with_order, `[[`, "order")
        I <- do.call(`[`, c(list(I), pos, list(drop = FALSE)))
    }

    .structure(list(domain = D,
                    incidence = I,
                    .arity = length(size),
                    .size = size),
               class = "relation_by_domain_and_incidence")
}

### ** .make_relation_by_domain_and_scores

.make_relation_by_domain_and_scores <-
function(D, scores)
{
    ## Assume a valid domain.

    n <- length(scores)

    .structure(list(domain = rep.int(list(D), 2L),
                    scores = scores,
                    .arity = 2L,
                    .size = c(n, n)),
               class = "relation_by_domain_and_scores")
}

### * Relation generators

### ** .make_relation_from_representation_and_meta

.make_relation_from_representation_and_meta <-
function(x, meta = NULL)
    .make_container(x, "relation", meta)

### ** .make_relation_from_domain_and_incidence

.make_relation_from_domain_and_incidence <-
function(D, I, meta = list())
{
    ## Canonicalize a valid incidence and determine whether it is crisp
    ## or not.
    I <- as.array(I)
    if(is.logical(I)) {
        I <- I + 0L
        is_crisp <- TRUE
    }
    else {
        ## Should be numeric if valid (maybe we should simply check this
        ## here as well?)
        is_crisp <- all((I %% 1) == 0)
    }
    R <- .make_relation_by_domain_and_incidence(D, I)
    meta["is_crisp"] <- is_crisp
    .make_relation_from_representation_and_meta(R, meta)
}

### ** .make_relation_from_domain_and_graph_components

.make_relation_from_domain_and_graph_components <-
function(D, G)
{
    ## (Assuming that G really has the graph *components* as obtained by
    ## .make_relation_graph_components().)
    values <- lapply(G, as.set)

    ## Get the domain.
    if(!is.null(D)) {
        L <- length(G)
        if((L > 0) && (length(D) != L))
            stop("Relation arity mismatch between domain and graph.")
        D <- lapply(D, as.cset)
        ## Check containment.
        ## <FIXME>
        ## Do we really need/want the factor transformation?
        ##    if((L > 0) && !all(mapply(set_is_subset,
        ##                              .transform_factors_into_characters(values),
        ##                              .transform_factors_into_characters(D))))
        ##         stop("Invalid graph with out-of-domain elements.")
        ## </FIXME>
        if((L > 0) && !all(mapply(set_is_subset, values, D)))
            stop("Invalid graph with out-of-domain elements.")
    } else
        D <- values

    ## Get the incidences.
    I  <- .make_incidence_from_domain_and_graph_components(D, G)
    M <- attr(G, "memberships")
    if (!is.null(M))
        I[I == 1L] <- M

    ## And put things together.
    .make_relation_from_domain_and_incidence(D, I)
}

### ** .make_relation_from_domain_and_scores

.make_relation_from_domain_and_scores <-
function(D, scores, meta = list())
{
    R <- .make_relation_by_domain_and_scores(D, scores)
    meta <- c(list(is_endorelation = TRUE, is_crisp = TRUE), meta)
    meta <- meta[!duplicated(names(meta))]
    .make_relation_from_representation_and_meta(R, meta)
}

### * Utilities

### ** .canonicalize_relation

## <NOTE>
## This should no longer be needed now that creating relations always
## canonicalizes to the unique set order.
## .canonicalize_relation <-
## function(R, D, pos = NULL)
## {
##     ## For a relation R with domain known to equal D in the sense that
##     ## the respective domain elements are the same sets (as tested for
##     ## by .domain_is_equal()), "canonicalize" R to have its domain
##     ## elements use the same *internal* order as the elements of D.
##     if(is.null(pos)) {
##         pos <- .match_domain_components(lapply(D, as.list), lapply(.domain(R), as.list))
##         ## If already canonical, do nothing.
##         if(!any(sapply(pos, is.unsorted))) return(R)
##     }
##     ## Use the reference domain, and reorder incidences.
##     I <- .reorder_incidence(.incidence(R), pos)
##     meta <- relation_properties(R)
##     .make_relation_from_domain_and_incidence(D, I, meta)
## }
## </NOTE>

### ** .make_data_frame_from_list

.make_data_frame_from_list <-
function(x, row.names = NULL)
{
    if(length(x)) {
        len <- length(x[[1L]])
        row.names <- if(!is.null(row.names)) {
            ## Do some checking ...
            if(length(row.names) != len)
                stop("Incorrect length of given 'row.names'.")
            as.character(row.names)
        }
        else
            .set_row_names(len)
    }
    .structure(x, class = "data.frame", row.names = row.names)
}

### * .make_domain_names_from_relation_graph_components

.make_domain_names_from_relation_graph_components <-
function(x, endorelation = FALSE) {
    ## In case the domains do not have names, use X_i.
    ## Needed for as.data.frame.relation();
    ## Not sure if we want this for relation_graph(), too.
    nms <- names(x)
    arity <- length(x)
    if (is.null(nms)) {
      nms <- if (endorelation) {
        rep.int("X", 2L)
        ## (Assuming that endorelations are always binary.)
      }
      else if (arity == 1L) "X"
      else sprintf("X_%d", seq_len(arity))
    }
    nms
}

### ** .make_incidence_from_domain_and_graph_components

.make_incidence_from_domain_and_graph_components <-
function(D, G, size = NULL)
{
    if(is.null(size)) size <- sapply(D, length)
    I <- array(0, size)
    if(length(G) > 0L)
        I[rbind(mapply(.exact_match,
                       G,
                       lapply(D, as.list)
                       )
                )
          ] <- 1
    I
}


### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***

