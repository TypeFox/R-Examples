### * relation_class_ids

relation_class_ids <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(!identical(relation_is_crisp(x), TRUE))
        stop("Argument 'x' must be a crisp relation with no missings.")
    if(relation_is_weak_order(x)) {
        ## Get the class ids of the corresponding indifference relation.
        ## One possibility:
        ##   I <- relation_incidence(x)
        ##   get_class_ids_from_incidence(I & t(I))
        ## But this is faster and also gives the class ids in the
        ## natural order of the indifference classes.
        s <- relation_scores(x, "ranks", decreasing = FALSE)
        ids <- match(s, sort(unique(s)))
        names(ids) <- names(s)
        ids
    }
    else if(relation_is_equivalence(x))
        get_class_ids_from_incidence(relation_incidence(x))
    else
        stop("Can only determine class ids for equivalences and weak orders.")
}

### * relation_classes

relation_classes <-
function(x)
{
    ids <- relation_class_ids(as.relation(x))
    out <- lapply(split(.get_elements_in_homorelation(x), ids), as.set)
    class(out) <- c("relation_classes_of_objects")
    out
}

print.relation_classes_of_objects <-
function(x, ...)
{
    y <- lapply(x, function(e) paste(format(e), collapse = "\n"))
    writeLines(formatDL(rev(names(x)), rev(unlist(y)),
                        style = "list", ...))
    invisible(x)
}


get_class_ids_from_incidence <-
function(x)
{
    ## Ugly ...
    y <- integer(nrow(x))
    c <- 1L
    pos <- seq_along(y)
    while(length(pos)) {
        ind <- x[pos[1L], pos] == 1
        y[pos[ind]] <- c
        pos <- pos[!ind]
        c <- c + 1L
    }
    names(y) <- rownames(x)
    y
}

### * relation_elements

## Let R be an endorelation.
## An element x is minimal if there is no "smaller" one, i.e.:
##   There is no y != x with y R x
## An element x is a first element if it is "not larger" than any
## other element, i.e.:
##   x R y for all y != x
## An element x is maximal if there is no "larger" one, i.e.:
##   There is no y != x with x R y
## An element x is a last element if it is "not smaller" than any
## other element, i.e.:
##   y R x for all y != x

## Note that sets cannot directly be indexed positionally.
## (Well, as of 2008-08-08 there is sets:::.set_subset() ...)

relation_elements <-
function(x, which, ...)
{
    ## We try to minimize code duplication, so that the respective
    ## getters only do computations based on incidences.  Of course, the
    ## 'elements' concept would at least make sense for arbitrary (not
    ## necessarily binary) homorelations ...
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Argument 'x' must be an endorelation.")
    I <- .incidence(x)
    if(any((I %% 1) != 0, na.rm = TRUE))
        stop("Argument 'x' must be a non-fuzzy relation.")
    
    which <- match.arg(which,
                       c("first", "minimal", "maximal", "last"))
    ## Allowing for user-defined 'which' functions is somewhat moot,
    ## because all we can do then is call which(x, ...), which users
    ## could have called directly in the first place.
    ind <- do.call(sprintf(".find_elements_being_%s", which),
                   list(I, ...))
    as.set(.get_elements_in_homorelation(x)[ind])
}

.find_elements_being_minimal <-
function(I, na.rm = TRUE)
{
    diag(I) <- 0
    colSums(I, na.rm = na.rm) == 0
}

.find_elements_being_first <-
function(I, na.rm = FALSE)
{
    diag(I) <- 1
    ind <- rowSums(I != 1, na.rm = na.rm) == 0
    !is.na(ind) & ind
}

.find_elements_being_last <-
function(I, na.rm = FALSE)
{
    diag(I) <- 1
    ind <- colSums(I != 1, na.rm = na.rm) == 0
    !is.na(ind) & ind
}

.find_elements_being_maximal <-
function(I, na.rm = TRUE)
{
    diag(I) <- 0
    rowSums(I, na.rm = na.rm) == 0
}

## * relation_successors

## See http://en.wikipedia.org/wiki/Covering_relation:
## Let X be a poset with associated strict partial order <.
## Then y covers x if x < y but there is no z with x < z and z < y.
## For the case of a general endorelation R, we use the following
## generalization: let P = R & !t(R) be the asymmetric part of R, and
## take y to cover x if x P y but there is no z with x P z and z P y.
##
## Note that there seem to be other definitions of covering as well.
## E.g., Brandt and Fischer, "Computational Aspects of Covering in
## Dominance Graphs" (In R. C. Holte and A. Howe, editors, Proceedings
## of the 22nd Conference on Artificial Intelligence (AAAI), pages
## 694-699. AAAI Press, 2007.), URL:
## http://www.tcs.informatik.uni-muenchen.de/~pamas/papers/aaai2007.pdf
## have upward, downward and bidirectional covering for dominance
## (asymmetric and irreflexive) relations > defined as
##   x C_u y: x > y and for all z, z > x implies z > y
##   x C_d y: x > y and for all z, y > z implies x > z
##   x C_b y: x C_u y and x C_d y
## Note that this is different from the Wikipedia definition.  If we
## take {1, 2, 3} with the natural > strict order and x=3 and y=1, then
## x C_b y because there is no z with z > x or y > z.
##
## Eventually, we could add a 'which' argument to relation_cover() ...

relation_successors <-
function(x, e = NULL)
{
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Argument 'x' must be an endorelation.")
    if(!identical(relation_is_crisp(x), TRUE))
        stop("Argument 'x' must be a crisp relation with no missings.")

    X <- .get_elements_in_homorelation(x)
    ## Need to find the positions of e in X.
    if(is.null(e)) {
        pos <- seq_along(X)
    } else {
        pos <- .exact_match(e, X)
        if(any(is.na(pos)))
            stop("Elements of 'e' must be contained in the domain components.")
    }
    ## Argh, terminology is really a nuisance.  If we had gone with the
    ## domain/codomain (domain/range) terminology commonly employed for
    ## endorelations, we could say that X is the domain of the relation,
    ## but then what is the tuple (X, X) called?  (And what is used in
    ## the general k-ary case?)  [Wikipedia says that X_1, ..., X_k are
    ## the domains of the relation.]
    I <- .relation_cover_incidences(relation_incidence(x), pos)
    out <- lapply(pos, function(p) as.set(X[I[p, ] == 1]))
    ## (Yes there may be more efficient ways ...)
    names(out) <- rownames(I)[pos]
    out
}

relation_precursors <-
function(x, e = NULL)
    relation_successors(t(x), e)

### * relation_cover

relation_cover <-
function(x)
{
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Argument 'x' must be an endorelation.")
    if(!identical(relation_is_crisp(x), TRUE))
        stop("Argument 'x' must be a crisp relation with no missings.")

    D <- .domain(x)
    I <- .relation_cover_incidences(.incidence(x), seq_along(D[[1L]]))
    meta <- list(is_endorelation = TRUE,
                 is_irreflexive = TRUE,
                 is_antisymmetric = TRUE,
                 is_asymmetric = TRUE)
    .make_relation_from_domain_and_incidence(D, I, meta)
}

.relation_cover_incidences <-
function(I, pos)
{
    ## Compute the incidences of the covering relation for the given
    ## positions.  When computing successors of a subset of elements
    ## this only computes what is needed (although not necessarily in
    ## the most effective way).

    ## Determine the incidences of the strict preference P(R) associated
    ## with R.
    I <- I * (1 - t(I))

    J <- I
    for(p in pos) {
        ## Compute all z for which x P(R) z.
        candidates <- which(I[p, ] == 1)
        ## Need to find those candidates y for which there is no z != y
        ## with z P(R) y.
        J[p, candidates] <-
            (colSums(I[candidates, candidates, drop = FALSE]) == 0)
    }
    J
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
