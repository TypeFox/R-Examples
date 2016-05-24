## Predicates.

## Note that we strongly prefer to use is_FOO rather than is.foo: we use
## this for distinguishing between class predicates (is.foo) and others.
## Alternatively, we could use a single predicate function
##   relation_test(x, predicate)
## The relation_is_${predicate}() approach has the advantage that some
## of these functions could be made generic eventually: provided we
## allow for relations without explicit coercion, we can simplify some
## of the tests (non-ordered factors given equivalence relations, etc).

### * Arity predicates

relation_is_binary <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    .arity(x) == 2L
}
relation_is_ternary <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    .arity(x) == 3L
}
relation_is_quaternary <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    .arity(x) == 4L
}

### * Predicates for general binary relations

## For now, always return FALSE if not crisp (with no NAs).

## Note that typically relations will have arity metadata, so binarity
## can be checked without computing incidences.

relation_is_left_total <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(!relation_is_binary(x)) return(FALSE)
    x <- relation_incidence(x)
    if(any(is.na(x)) || !all((x %% 1) == 0)) return(FALSE)
    all(rowSums(x) >= 1)
}

relation_is_right_total <-
relation_is_surjective <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(!relation_is_binary(x)) return(FALSE)
    x <- relation_incidence(x)
    if(any(is.na(x)) || !all((x %% 1) == 0)) return(FALSE)
    all(colSums(x) >= 1)
}

relation_is_functional <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(!relation_is_binary(x)) return(FALSE)
    x <- relation_incidence(x)
    if(any(is.na(x)) || !all((x %% 1) == 0)) return(FALSE)
    all(rowSums(x) <= 1)
}

relation_is_injective <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(!relation_is_binary(x)) return(FALSE)
    x <- relation_incidence(x)
    if(any(is.na(x)) || !all((x %% 1) == 0)) return(FALSE)
    all(colSums(x) <= 1)
}

relation_is_bijective <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(!relation_is_binary(x)) return(FALSE)
    x <- relation_incidence(x)
    if(any(is.na(x)) || !all((x %% 1) == 0)) return(FALSE)
    all(rowSums(x) == 1) && all(colSums(x) == 1)
}

### * Endorelations and predicates of such

relation_is_endorelation <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(.has_property(x, "is_endorelation"))
        return(.get_property(x, "is_endorelation"))
    .relation_is_endorelation_using_incidence(relation_incidence(x))
}

.relation_is_endorelation_using_incidence <-
function(x)
{
    ## For internal purposes only: used to avoid the need for possibly
    ## computing incidences at least twice in some of the code below
    ## (not that much of an issue as long as incidences are the only
    ## possible representation).

    ## Need some heuristic to determine whether we have an endorelation
    ## or not.  Idea: assume yes if nrow = ncol and either there are no
    ## dimnames (argh) or rownames and colnames are identical (better);
    ## otherwise, assume no.
    (is.matrix(x)
     && (nrow(x) == ncol(x))
     && identical(rownames(x), colnames(x)))
}

relation_is_homogeneous <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(.has_property(x, "is_homogeneous"))
        return(.get_property(x, "is_homogeneous"))
    length(unique(.domain(x))) == 1L
}

## <NOTE>
## When inferring predicates from properties, assume that these are
## computed correctly for the default case (na.rm = FALSE).  Otherwise,
## note that predicates are always computed from conjunctive normal
## forms
##   all(c(test_1, ..., test_k))
## where na.rm = TRUE removes the NA tests.  Hence, NA removal can never
## change FALSE to TRUE.
## </NOTE>

relation_is_crisp <-
function(x, na.rm = FALSE)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(.has_property(x, "is_crisp")) {
        p <- .get_property(x, "is_crisp")
        if(!na.rm || identical(p, FALSE))
            return(p)
    }
    all((relation_incidence(x) %% 1) == 0, na.rm = na.rm)
}

## <NOTE>
## Sometimes "total" is used synonymously to "complete".
## http://en.wikipedia.org/wiki/Binary_relation has two different usages
## of "total" ...
## </NOTE>
relation_is_complete <-
function(x, na.rm = FALSE)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(.has_property(x, "is_complete")) {
        p <- .get_property(x, "is_complete")
        if(!na.rm || identical(p, FALSE))
            return(p)
    }
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    ## Use
    ##   all(.S.(x, t(x))[row(x) != col(x)] == 1, na.rm = na.rm)
    ## or more efficiently:
    diag(x) <- 1
    all(.S.(x, t(x)) == 1, na.rm = na.rm)
}

relation_is_match <-
relation_is_strongly_complete <-
function(x, na.rm = FALSE)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    ## We currently never set is_strongly_complete metadata.
    ## We could look at both is_complete and is_reflexive, though.
    all(.S.(x, t(x)) == 1, na.rm = na.rm)
}

relation_is_reflexive <-
function(x, na.rm = FALSE)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(.has_property(x, "is_reflexive")) {
        p <- .get_property(x, "is_reflexive")
        if(!na.rm || identical(p, FALSE))
            return(p)
    }
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    all(diag(x) == 1, na.rm = na.rm)
}

relation_is_irreflexive <-
function(x, na.rm = FALSE)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    all(diag(x) == 0, na.rm = na.rm)
}

relation_is_coreflexive <-
function(x, na.rm = FALSE)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    ## Use
    ##   all(x[row(x) != col(x)] == 0, na.rm = na.rm))
    ## or more efficiently:
    diag(x) <- 0
    all(x == 0, na.rm = na.rm)
}

relation_is_symmetric <-
function(x, na.rm = FALSE)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(.has_property(x, "is_symmetric")) {
        p <- .get_property(x, "is_symmetric")
        if(!na.rm || identical(p, FALSE))
            return(p)
    }
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    ## Use
    ##   all((x == t(x))[row(x) != col(x)], na.rm = na.rm)
    ## or more efficiently:
    diag(x) <- 1
    all(x == t(x), na.rm = na.rm)
}

relation_is_asymmetric <-
function(x, na.rm = FALSE)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    all(.T.(x, t(x)) == 0, na.rm = na.rm)
}

relation_is_antisymmetric <-
function(x, na.rm = FALSE)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(.has_property(x, "is_antisymmetric")) {
        p <- .get_property(x, "is_antisymmetric")
        if(!na.rm || identical(p, FALSE))
            return(p)
    }
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    ## Use
    ##   all(.T.(x, t(x))[row(x) != col(x)] == 0, na.rm = na.rm)
    ## or more efficiently:
    diag(x) <- 0
    all(.T.(x, t(x)) == 0, na.rm = na.rm)
}

relation_is_transitive <-
function(x, na.rm = FALSE)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(.has_property(x, "is_transitive")) {
        p <- .get_property(x, "is_transitive")
        if(!na.rm || identical(p, FALSE))
            return(p)
    }
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    p <- TRUE
    if(na.rm) {
        for(j in seq_len(ncol(x)))
            if(any(outer(x[, j], x[j, ], .T.) > x, na.rm = TRUE))
                return(FALSE)
    } else {
        for(j in seq_len(ncol(x))) {
            v <- any(outer(x[, j], x[j, ], .T.) > x)
            if(is.na(v)) p <- NA
            else if(v) return(FALSE)
        }
    }
    p
}

relation_is_negatively_transitive <-
function(x, na.rm = FALSE)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    p <- TRUE
    if(na.rm) {
        for(j in seq_len(ncol(x)))
            if(any(outer(x[, j], x[j, ], .S.) < x, na.rm = TRUE))
                return(FALSE)
    } else {
        for(j in seq_len(ncol(x))) {
            v <- any(outer(x[, j], x[j, ], .S.) < x)
            if(is.na(v)) p <- NA
            else if(v) return(FALSE)
        }
    }
    p
}

relation_is_quasitransitive <-
function(x, na.rm = FALSE)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    ## In essence, once we know that x is an endorelation, we can do
    ##   relation_is_transitive(x & !t(x))
    ## But to do so, we need to look at the incidences ...
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    x <- .T.(x, .N.(t(x)))
    p <- TRUE
    if(na.rm) {
        for(j in seq_len(ncol(x)))
            if(any(outer(x[, j], x[j, ], .T.) > x, na.rm = TRUE))
                return(FALSE)
    } else {
        for(j in seq_len(ncol(x))) {
            v <- any(outer(x[, j], x[j, ], .T.) > x)
            if(is.na(v)) p <- NA
            else if(v) return(FALSE)
        }
    }
    p
}

relation_is_Ferrers <-
function(x, na.rm = FALSE)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    n <- nrow(x)
    p <- TRUE
    ## Relation is Ferrers iff for all i, j, k, l
    ##   T(R(i,j), R(k,l)) <= S(R(i,l), R(k,j))
    if(na.rm) {
        for(j in seq_len(n))
            for(l in seq_len(n)) {
                if(any(outer(x[, j], x[, l], .T.) >
                       outer(x[, l], x[, j], .S.),
                       na.rm = TRUE))
                    return(FALSE)
            }
    } else {
        for(j in seq_len(n))
            for(l in seq_len(n)) {
                v <- any(outer(x[, j], x[, l], .T.) >
                         outer(x[, l], x[, j], .S.))
                if(is.na(v)) p <- NA
                else if(v) return(FALSE)
            }
    }
    p
}

relation_is_semitransitive <-
function(x, na.rm = FALSE)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    n <- nrow(x)
    p <- TRUE
    ## Relation is semitransitive iff for all i, j, k, l
    ##   T(R(i,k), R(k,j)) <= S(R(i,l), R(l,j))
    if(na.rm) {
        for(k in seq_len(n))
            for(l in seq_len(n)) {
                if(any(outer(x[, k], x[k, ], .T.) >
                       outer(x[, l], x[l, ], .S.),
                       na.rm = TRUE))
                    return(FALSE)
            }
    } else {
        for(k in seq_len(n))
            for(l in seq_len(n)) {
                v <- any(outer(x[, k], x[k, ], .T.) >
                         outer(x[, l], x[l, ], .S.))
                if(is.na(v)) p <- NA
                else if(v) return(FALSE)
            }
    }
    p
}

relation_is_trichotomous <-
function(x, na.rm = FALSE)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    ## Trichotomous: exactly one of xRy, yRx, or x=y holds.
    ## Interpret this literally so that for y=x xRx cannot hold (i.e.,
    ## the relation must be irreflexive).
    (all(diag(x) == 0, na.rm = na.rm) &&
     all(abs(x - t(x))[row(x) != col(x)] == 1, na.rm = na.rm))
}

relation_is_Euclidean <-
function(x, na.rm = FALSE)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    p <- TRUE
    ## Euclidean: xRy & xRz => yRz.
    if(na.rm) {
        for(j in seq_len(nrow(x)))
            if(any(outer(x[j, ], x[j, ], .T.) > x, na.rm = TRUE))
                return(FALSE)
    } else {
        for(j in seq_len(nrow(x))) {
            v <- any(outer(x[j, ], x[j, ], .T.) > x)
            if(is.na(v)) p <- NA
            else if(v) return(FALSE)
        }
    }
    p
}

## And now combine:
## Of course, these could be made more efficient by doing all
## computations [on incidences] just once ...

relation_is_equivalence <-
function(x, na.rm = FALSE)
    (relation_is_endorelation(x)
     && relation_is_reflexive(x, na.rm)
     && relation_is_symmetric(x, na.rm)
     && relation_is_transitive(x, na.rm))

relation_is_weak_order <-
relation_is_preference <-
function(x, na.rm = FALSE)
{
    ## Equivalently: strongly complete and transitive.
    (relation_is_endorelation(x)
     && relation_is_complete(x, na.rm)
     && relation_is_reflexive(x, na.rm)
     && relation_is_transitive(x, na.rm))
}

relation_is_preorder <-
relation_is_quasiorder <-
function(x, na.rm = FALSE)
    (relation_is_endorelation(x)
     && relation_is_reflexive(x, na.rm)
     && relation_is_transitive(x, na.rm))

relation_is_partial_order <-
function(x, na.rm = FALSE)
    (relation_is_endorelation(x)
     && relation_is_reflexive(x, na.rm)
     && relation_is_antisymmetric(x, na.rm)
     && relation_is_transitive(x, na.rm))

relation_is_linear_order <-
function(x, na.rm = FALSE)
    (relation_is_partial_order(x, na.rm)
     && relation_is_complete(x, na.rm))

relation_is_strict_partial_order <-
function(x, na.rm = FALSE)
    (relation_is_endorelation(x)
     && relation_is_asymmetric(x, na.rm)
     && relation_is_transitive(x, na.rm))

relation_is_strict_linear_order <-
function(x, na.rm = FALSE)
    (relation_is_strict_partial_order(x, na.rm)
     && relation_is_complete(x, na.rm))

relation_is_tournament <-
function(x, na.rm = FALSE)
{
    ## The references disagree on whether tournaments are reflexive
    ## (e.g., Barthelemy) or irreflexive (e.g., Fodor/Roubens).
    ## We use the latter (as implied by asymmetry).
    (relation_is_endorelation(x)
     && relation_is_complete(x, na.rm)
     && relation_is_asymmetric(x, na.rm))
}

relation_is_interval_order <-
function(x, na.rm = FALSE)
    (relation_is_endorelation(x)
     && relation_is_complete(x, na.rm)
     && relation_is_Ferrers(x, na.rm))

relation_is_semiorder <-
function(x, na.rm = FALSE)
    (relation_is_interval_order(x, na.rm)
     && relation_is_semitransitive(x, na.rm))

relation_has_missings <-
function(x)
    any(is.na(.incidence(x)))

relation_is_acyclic <-
function(x)
    relation_is_antisymmetric(transitive_closure(x))

relation_is_cyclic <- 
function(x)
    !relation_is_acyclic(x)

## FIXME: add meta cache to improve performance
## (idea: all predicates should call .foo-methods internally with additional
## meta-argument, where all intermediate results are stored and which is
## looked up before computations.)
.check_all_predicates <-
function(x, ...)
{
    preds <- ls("package:relations", pattern = "relation_is_.*")
    props <- sapply(preds, do.call, c(list(x), list(...)))
    names(props) <- sub("relation_is_", "", names(props))
    props
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
