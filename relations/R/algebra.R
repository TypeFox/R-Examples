### Relational algebra-like operations.

### * relation_projection

relation_projection <-
function(x, margin = NULL)
{
    if (is.null(margin))
        return(x)
    D <- relation_domain(x)

    ind <- if (is.character(margin))
        match(margin, names(D))
    else if (is.numeric(margin))
        match(margin, seq_along(D))
    else NULL

    if (!length(ind) || any(is.na(ind)))
        stop("Invalid projection margin.")

    ## Generally, incidences are aggregated using the T-conorm.
    ## However, for crisp relations, we use 'any' for performance reasons.
    I <- relation_incidence(x)
    S.FUN <- if(isTRUE(relation_is_crisp(x))) {
        mode(I) <- "logical"
        any
    } else
        function(i) Reduce(.S., i)
    .make_relation_from_domain_and_incidence(.domain(x)[ind],
                                             apply(I, ind, S.FUN))
}

aperm.relation <-
function(a, perm = NULL, ...)
{
    D <- relation_domain(a)
    s <- seq_along(D)

    if(is.null(perm)) 
        ind <- rev(s)
    else {
        ind <- if(is.character(perm))
            match(perm, names(D))
        else if(is.numeric(perm))
            match(perm, s)
        else
            NULL
        if(length(ind) != length(s) || any(is.na(ind)))
            stop("Invalid permutation.")
    }

    I <- relation_incidence(a)
    .make_relation_from_domain_and_incidence(.domain(a)[ind],
                                             aperm(I, ind))
}

### * relation_selection (= subset)

## <FIXME>
## Do we really have to copy the code from subset.data.frame()?
## </FIXME>
relation_selection <-
function(x, subset)
{
    .stop_if_not_relation_has_unique_domain_names(x)
    df <- as.data.frame(x)
    if (missing(subset))
        r <- TRUE
    else {
        e <- substitute(subset)
        r <- eval(e, df, parent.frame())
        if (!is.logical(r))
            stop("'subset' must evaluate to logical")
        r <- r & !is.na(r)
    }
    relation_graph(x) <- df[r,]
    x
}

"%><%" <-
relation_cartesian <-
function(x, y, ...)
{
    if (missing(y)) return(x)
    l <- list(...)
    if(length(l))
        return(Recall(x, Recall(y, ...)))
    T.FUN <- if(relation_is_crisp(x, na.rm = TRUE) &&
                relation_is_crisp(y, na.rm = TRUE))
        "*"
    else
        .T.
    .make_relation_from_domain_and_incidence(c(.domain(x), .domain(y)),
                                             outer(.incidence(x),
                                                   .incidence(y),
                                                   T.FUN)
                                             )
}

### * relation_union

## For union of relations, we also allow non-identical domain labels.

"%U%" <-
relation_union <-
function(x, y, ...)
{
    ## handle 1 and more than 2 arguments
    if (missing(y)) return(x)
    l <- list(...)
    if(length(l))
        return(Recall(x, Recall(y, ...)))

    ## check arities
    ## (use relation_domain() instead of .domain() here,
    ## since we need tuples of sets for the combination.)
    Dx <- relation_domain(x)
    Dy <- relation_domain(y)
    if (length(Dx) != length(Dy))
        stop("Relation arity mismatch.")

    ## combine domains
    Dxy <- Map(set_union, Dx, Dy)

    ## extract incidences for combined domain
    Ix <- Iy <- array(0, sapply(Dxy, length),
                      lapply(Dxy, LABELS, quote = FALSE))
    Ix <- do.call("[<-", c(list(Ix),
                           lapply(Dx, LABELS, quote = FALSE),
                           list(relation_incidence(x))))
    Iy <- do.call("[<-", c(list(Iy),
                           lapply(Dy, LABELS, quote = FALSE),
                           list(relation_incidence(y))))

    ## and put things together
    .make_relation_from_domain_and_incidence(Dxy, .S.(Ix, Iy))
}


### * relation_intersection

## For the intersection, we also allow non-identical domain levels.

relation_intersection <-
function(x, y, ...)
{
    ## handle 1 and more than 2 arguments
    if (missing(y)) return(x)
    l <- list(...)
    if(length(l))
        return(Recall(x, Recall(y, ...)))

    ## check arities
    ## (use relation_domain() instead of .domain() here,
    ## since we need tuples of sets for the combination.)
    Dx <- relation_domain(x)
    Dy <- relation_domain(y)
    if (length(Dx) != length(Dy))
        stop("Relation arity mismatch.")

    ## intersect domains
    Dxy <- Map(set_intersection, Dx, Dy)

    ## check if there is some overlap
    if (any(sapply(Dxy, set_is_empty)))
        return(set())

    ## extract incidences for common domain
    Ix <- do.call("[", c(list(relation_incidence(x)),
                         lapply(Dxy, LABELS, quote = FALSE)))
    Iy <- do.call("[", c(list(relation_incidence(y)),
                         lapply(Dxy, LABELS, quote = FALSE)))

    ## and put things together
    .make_relation_from_domain_and_incidence(Dxy, .T.(Ix, Iy))
}

### * relation_complement

relation_complement <-
function(x, y = NULL)
{
    ## handle unary case
    if (is.null(y))
        return(.make_relation_from_domain_and_incidence(.domain(x),
                                                        .N.(.incidence(x))))

    Dx <- .domain(x)
    Dy <- .domain(y)
    if (length(Dx) != length(Dy))
        stop("Relation arity mismatch.")

    ## extract incidences for common domain
    I <- do.call("[", c(list(.incidence(y)), Map(.exact_match, Dx, Dy)))
    I[is.na(I)] <- 0
    relation_intersection(x,
        relation_complement(.make_relation_from_domain_and_incidence(Dx, I)))
}

### * relation_symdiff

relation_symdiff <-
function(x, y)
    relation_union(relation_complement(x, y),
                   relation_complement(y, x))

### * relation_division

relation_division <-
function(x, y)
{
    .stop_if_not_relation_has_unique_domain_names(x)
    .stop_if_not_relation_has_unique_domain_names(y)

    if (length(relation_graph(y)) < 1L)
        stop("Division by empty relations not defined.")

    dx <- relation_domain_names(x)
    dy <- relation_domain_names(y)

    if (!all(dy %in% dx))
        stop("Divisor domain must be a subset of the dividend domain.")

    ## find attributes unique to x
    dxunique <- dx[!dx %in% dy]
    if (length(dxunique) < 1L)
        stop("Dividend needs at least one unique domain.")

    ## create projection of x to its unique attributes
    xunique <- relation_projection(x, dxunique)

    ## compute "maximum" set of tuples
    T <- relation_cartesian(xunique, y)

    ## remove actual set of tuples, and remove the projection
    ## of the remaining sets from the dividend
    relation_complement(xunique,
                        relation_projection(relation_complement(T, x),
                                            dxunique))
}

### relation_remainder

relation_remainder <-
function(x, y)
    relation_complement(x,
                        relation_cartesian(relation_division(x, y), y))

### * relation_join et al

"%|><|%" <-
relation_join <-
function(x, y, ...)
{
    ## check domains
    .stop_if_not_relation_has_unique_domain_names(x)
    .stop_if_not_relation_has_unique_domain_names(y)

    ## add memberships, if any
    X <- as.data.frame(x)
    Y <- as.data.frame(y)
    nms <- unique(c(names(X), names(Y)))
    fuzzy <- !isTRUE(relation_is_crisp(x)) || !isTRUE(relation_is_crisp(y))
    if (fuzzy) {
        Mx <- attr(X, "memberships")
        if (is.null(Mx)) Mx <- 1
        X <- cbind(X, "MEMBERSHIP.x" = Mx)

        My <- attr(Y, "memberships")
        if (is.null(My)) My <- 1
        Y <- cbind(Y, "MEMBERSHIP.y" = My)
    }

    ## use merge for the operation
    tmp <- merge(X, Y, ...)

    ## handle empty result
    if(nrow(tmp) < 1L)
        return(set())

    ## combine memberships for fuzzy relations
    M <- if (fuzzy) {
        Mx <- tmp[,"MEMBERSHIP.x"]
        Mx[is.na(Mx)] <- 1
        My <- tmp[,"MEMBERSHIP.y"]
        My[is.na(My)] <- 1
        .T.(Mx, My)
    } else
        NULL

    ## rearrange columns & return relation
    as.relation(.structure(tmp[, nms], memberships = M))
}

"%><=%" <-
function(x, y, ...)
    relation_join(x, y, all.y = TRUE, ...)

"%=><%" <-
function(x, y, ...)
    relation_join(x, y, all.x = TRUE, ...)

"%=><=%" <-
function(x, y, ...)
    relation_join(x, y, all = TRUE, ...)

"%|><%" <-
relation_semijoin <-
function(x, y, ...)
    relation_projection(relation_join(x, y, ...),
                        relation_domain_names(x))

"%><|%" <-
function(x, y, ...)
    relation_semijoin(y, x, ...)

"%|>%" <-
relation_antijoin <-
function(x, y, ...)
    x - relation_semijoin(x, y, ...)

"%<|%" <-
function(x, y, ...)
    relation_antijoin(y, x, ...)

### * .stop_if_not_relation_has_unique_domain_names

.stop_if_not_relation_has_unique_domain_names <-
function(x)
{
    nms <- relation_domain_names(x)
    if(is.null(nms) || (length(nms) < .arity(x)) || any(duplicated(nms)))
        stop("Relation(s) with unique domain names required.")
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
