relation_violations <-
function(x,
         property =
         c("complete", "match",
           "reflexive", "irreflexive", "coreflexive",
           "symmetric", "antisymmetric", "asymmetric",
           "transitive", "negatively_transitive", "Ferrers",
           "semitransitive",
           "trichotomous",
           "Euclidean"),
         tuples = FALSE,
         na.rm = FALSE)
{
    if (!relation_is_endorelation(x))
        stop("Relation violations only defined for endorelations.")

    property <- match.arg(property)
    I <- .incidence(x)

    if(!tuples) {
        do.call(sprintf(".amount_by_which_relation_is_not_%s", property),
                list(I, na.rm = na.rm))
    } else {
        ## First get a matrix of indices of violating tuples.
        ind <- do.call(sprintf(".tuples_for_which_relation_is_not_%s",
                               property),
                       list(I, na.rm = na.rm))
        if(!nrow(ind)) return(set())
        ## And construct a set of violating tuples from this.
        D <- rep.int(list(.get_elements_in_homorelation(x)), ncol(ind))
        as.set(apply(ind, 1L, function(i) as.tuple(Map(`[[`, D, i))))
    }
}

.amount_by_which_relation_is_not_complete <-
function(I, na.rm = FALSE)
{
    diag(I) <- 1
    sum(1 - .S.(I, t(I)), na.rm = na.rm) / 2
}

.tuples_for_which_relation_is_not_complete <-
function(I, na.rm = FALSE)
{
    diag(I) <- 1
    ind <- .which(.S.(I, t(I)) < 1, arr.ind = TRUE, na.rm = na.rm)
    ind[ind[, 1L] <= ind[, 2L], , drop = FALSE]
}

.amount_by_which_relation_is_not_match <-
function(I, na.rm = FALSE)
{
    I <- 1 - .S.(I, t(I))
    D <- diag(I)
    diag(I) <- 0
    sum(I, na.rm = na.rm) / 2 + sum(D, na.rm = na.rm)
}

.tuples_for_which_relation_is_not_match <-
function(I, na.rm = FALSE)
{
    I <- .S.(I, t(I))
    ind <- .which(diag(I) < 1, na.rm = na.rm)
    diag(I) <- 1
    ind <- rbind(cbind(ind, ind),
                 .which(I < 1, arr.ind = TRUE, na.rm = na.rm))
    ind[ind[, 1L] <= ind[, 2L], , drop = FALSE]
}

.amount_by_which_relation_is_not_reflexive <-
function(I, na.rm = FALSE)
    sum(1 - diag(I), na.rm = na.rm)

.tuples_for_which_relation_is_not_reflexive <-
function(I, na.rm = FALSE)
    matrix(.which(diag(I) < 1, na.rm = na.rm), ncol = 1L)

.amount_by_which_relation_is_not_irreflexive <-
function(I, na.rm = FALSE)
    sum(diag(I), na.rm = na.rm)

.tuples_for_which_relation_is_not_irreflexive <-
function(I, na.rm = FALSE)
    matrix(.which(diag(I) > 0, na.rm = na.rm), ncol = 1L)
    
.amount_by_which_relation_is_not_coreflexive <-
function(I, na.rm = FALSE)
{
    diag(I) <- 0
    sum(I, na.rm = na.rm)
}

.tuples_for_which_relation_is_not_coreflexive <-
function(I, na.rm = FALSE)
{
    diag(I) <- 0
    .which(I > 0, arr.ind = TRUE, na.rm = na.rm)
}

.amount_by_which_relation_is_not_symmetric <-
function(I, na.rm = FALSE)
{
    diag(I) <- 1
    sum(abs(I - t(I)), na.rm = na.rm) / 2
}

.tuples_for_which_relation_is_not_symmetric <-
function(I, na.rm = FALSE)
{
    diag(I) <- 1
    .which(I != t(I), arr.ind = TRUE, na.rm = na.rm)
}

.amount_by_which_relation_is_not_asymmetric <-
function(I, na.rm = FALSE)
{
    I <- .T.(I, t(I))
    D <- diag(I)
    diag(I) <- 0
    sum(I, na.rm = na.rm) / 2 + sum(D, na.rm = na.rm)
}

.tuples_for_which_relation_is_not_asymmetric <-
function(I, na.rm = FALSE)
{
    ind <- .which(.T.(I, t(I)) > 0, arr.ind = TRUE, na.rm = na.rm)
    ind[ind[, 1L] <= ind[, 2L], , drop = FALSE]
}

.amount_by_which_relation_is_not_antisymmetric <-
function(I, na.rm = FALSE)
{
    diag(I) <- 0
    sum(.T.(I, t(I)), na.rm = na.rm) / 2
}

.tuples_for_which_relation_is_not_antisymmetric <-
function(I, na.rm = FALSE)
{
    ind <- .which(.T.(I, t(I)) > 0, arr.ind = TRUE, na.rm = na.rm)
    ind[ind[, 1L] < ind[, 2L], , drop = FALSE]
}    

.amount_by_which_relation_is_not_transitive <-
function(I, na.rm = FALSE)
    sum(sapply(seq_len(nrow(I)),
               function(j) pmax(outer(I[, j], I[j, ], .T.) - I, 0,
                                na.rm = na.rm)))

.tuples_for_which_relation_is_not_transitive <-
function(I, na.rm = FALSE)
{
    ind <- lapply(seq_len(nrow(I)),
                  function(j) {
                      pos <- .which(outer(I[, j], I[j, ], .T.) > I,
                                    arr.ind = TRUE, na.rm = na.rm)
                      cbind(pos, rep.int(j, nrow(pos)))
                  })
    do.call(rbind, ind)[, c(1L, 3L, 2L), drop = FALSE]
}

.amount_by_which_relation_is_not_negatively_transitive <-
function(I, na.rm = FALSE)
    sum(sapply(seq_len(nrow(I)),
               function(j) pmax(I - outer(I[, j], I[j, ], .S.), 0,
                                na.rm = na.rm)))

.tuples_for_which_relation_is_not_negatively_transitive <-
function(I, na.rm = FALSE)
{
    ind <- lapply(seq_len(nrow(I)),
                  function(j) {
                      pos <- .which(outer(I[, j], I[j, ], .S.) < I,
                                    arr.ind = TRUE, na.rm = na.rm)
                      cbind(pos, rep.int(j, nrow(pos)))
                  })
    do.call(rbind, ind)[, c(1L, 3L, 2L), drop = FALSE]
}

.amount_by_which_relation_is_not_Ferrers <-
function(I, na.rm = FALSE)
{
    out <- 0
    for(j in seq_len(nrow(I)))
        for(l in seq_len(nrow(I)))
            out <- out + sum(pmax(outer(I[, j], I[, l], .T.) -
                                  outer(I[, l], I[, j], .S.),
                                  0,
                                  na.rm = na.rm))
    out
}

.tuples_for_which_relation_is_not_Ferrers <-
function(I, na.rm = FALSE)
{
    n <- nrow(I)
    ind <- Map(function(j, l) {
                   pos <- .which(outer(I[, j], I[, l], .T.) >
                                 outer(I[, l], I[, j], .S.),
                                 arr.ind = TRUE, na.rm = na.rm)
                   cbind(pos,
                         rep.int(j, nrow(pos)),
                         rep.int(l, nrow(pos)))
               },
               rep.int(seq_len(n), n),
               rep(seq_len(n), each = n))
    do.call(rbind, ind)[, c(1L, 3L, 2L, 4L), drop = FALSE]
}

.amount_by_which_relation_is_not_semitransitive <-
function(I, na.rm = FALSE)
{
    out <- 0
    for(k in seq_len(nrow(I)))
        for(l in seq_len(nrow(I)))
            out <- out + sum(pmax(outer(I[, k], I[k, ], .T.) -
                                  outer(I[, l], I[l, ], .S.),
                                  0,
                                  na.rm = na.rm))
    out
}

.tuples_for_which_relation_is_not_semitransitive <-
function(I, na.rm = FALSE)
{
    n <- nrow(I)
    ind <- Map(function(k, l) {
                   pos <- .which(outer(I[, k], I[k, ], .T.) >
                                 outer(I[, l], I[l, ], .S.),
                                 arr.ind = TRUE, na.rm = na.rm)
                   cbind(pos,
                         rep.int(k, nrow(pos)),
                         rep.int(l, nrow(pos)))
               },
               rep.int(seq_len(n), n),
               rep(seq_len(n), each = n))
    do.call(rbind, ind)
}   

.amount_by_which_relation_is_not_trichotomous <-
function(I, na.rm = FALSE)
    (sum(diag(I), na.rm = na.rm) +
     sum(1 - abs(I - t(I))[row(I) != col(I)], na.rm = na.rm) / 2)

.tuples_for_which_relation_is_not_trichotomous <-
function(I, na.rm = FALSE)
{
    ind <- .which(abs(I - t(I)) < 1, arr.ind = TRUE, na.rm = na.rm)
    ind <- ind[ind[, 1L] < ind[, 2L], , drop = FALSE]
    pos <- .which(diag(I) > 0, na.rm = na.rm)
    rbind(ind, cbind(pos, pos))
}

.amount_by_which_relation_is_not_Euclidean <-
function(I, na.rm = FALSE)
    sum(sapply(seq_len(nrow(I)),
               function(i) pmax(outer(I[i, ], I[i, ], .T.) - I, 0,
                                na.rm = na.rm)))

.tuples_for_which_relation_is_not_Euclidean <-
function(I, na.rm = FALSE)
{
    ind <- lapply(seq_len(nrow(I)),
                  function(i) {
                      pos <- .which(outer(I[i, ], I[i, ], .T.) > I,
                                    arr.ind = TRUE, na.rm = na.rm)
                      cbind(rep.int(i, nrow(pos)), pos)
                  })
    do.call(rbind, ind)
}

.which <-
function(x, arr.ind = FALSE, na.rm = TRUE)
{
    if(!na.rm) x <- x | is.na(x)
    which(x, arr.ind = arr.ind)
}
