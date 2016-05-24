### Consensus relations.

### * relation_consensus

relation_consensus <-
function(x, method = NULL, weights = 1, control = list(), ...)
{
    dots <- list(...)
    control[names(dots)] <- dots

    if(inherits(x, "gset")) {
        relations <- relation_ensemble(list = gset_support(x))
        weights <- gset_memberships(x)
    } else {
        relations <- as.relation_ensemble(x)
    }

    if(!length(relations))
        stop("Cannot compute consensus of empty ensemble.")

    weights <- rep(weights, length.out = length(relations))
    if(any(weights < 0))
        stop("Argument 'weights' has negative elements.")
    if(!any(weights > 0))
        stop("Argument 'weights' has no positive elements.")

    if(!is.function(method)) {
        if(!inherits(method, "relation_consensus_method")) {
            ## Get the method definition from the registry.
            if(!is.character(method) || (length(method) != 1L))
                stop("Invalid 'method' argument.")
            entry <- get_relation_consensus_method(method)
            if(is.null(entry))
                stop(gettextf("Method '%s' is not a valid consensus method.",
                              method),
                     domain = NA)
            method <- entry
        }
        method <- method$definition
    }

    method(relations, weights, control)
}

### * Relation consensus "methods"

### ** .relation_consensus_Borda

## Kendall/Borda method
.relation_consensus_Borda <-
function(relations, weights, control)
    .relation_consensus_score(relations, weights, control,
                              relation_scores, "Borda")

### ** .relation_consensus_Copeland

## Copeland method
.relation_consensus_Copeland <-
function(relations, weights, control)
    .relation_consensus_score(relations, weights, control,
                              relation_scores, "Copeland")

### ** .relation_consensus_score

.relation_consensus_score <-
function(relations, weights, control, FUN, ...)
{
    ## Several sanity checks could be done here.
    ## In particular, check whether all relations are in fact complete
    ## preferences (or whatever is really necessary).

    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")

    ## extract options.
    n <- .n_from_control_list(control)
    f_l_o <- isTRUE(control$L) ## force linear order?

    ## Get the scores.
    scores <- lapply(relations, FUN, ...)

    ## Multiply by the weights and compute the total scores.
    scores <- rowSums(mapply("*", scores, weights))

    ## break ties directly if a single linear order is enforced.
    out <-
        rank(scores,
             ties.method = if (f_l_o && n == 1L) "first" else "average")

    INC <- function(S) outer(S, S, "<=")
    I <- if (f_l_o && n > 1L) {
        ## find all groupwise permutations & combine
        l <- expand.grid(lapply(split(seq_along(out), out), .permute))

        ## Only use up to n combinations
        l <- l[seq_len(min(n, nrow(l))),,drop = FALSE]

        ## create incidences
        unlist(apply(l, 1, function(i) list(INC(unlist(i)))), recursive = FALSE)
    } else INC(out)

    meta <- list(is_endorelation = TRUE,
                 is_complete = TRUE,
                 is_reflexive = TRUE,
                 is_antisymmetric = f_l_o || !any(duplicated(out)),
                 is_transitive = TRUE,
                 scores = scores)
    .make_consensus_from_incidences(.domain(relations), I, meta)
}

### ** .relation_consensus_majority

.relation_consensus_majority <-
function(relations, weights, control)
{
    p <- control$p
    if(is.null(p)) p <- 1 / 2
    ## Add some sanity checking for p eventually.

    incidences <- lapply(relations, relation_incidence)
    weights <- rep(weights, length.out = length(relations))

    I <- (.weighted_sum_of_arrays(incidences, weights, na.rm = TRUE) /
          .weighted_sum_of_arrays(lapply(incidences, is.finite),
                                  weights))
    I <- if(p == 1)
        I == 1
    else
        I > p

    ## (Could add a tie-breaking mechanism a la Condorcet, with an
    ## option for finding all solutions.)

    .make_relation_from_domain_and_incidence(.domain(relations), I)
}

### ** .relation_consensus_CS

.relation_consensus_CS <-
function(relations, weights, control)
{
    ## Cook and Seiford, Management Science (1978).
    ## Determine a linear order minimizing the aggregate Cook-Seiford
    ## dissimilarity to a given ensemble of relations (originally:
    ## complete rankings [i.e., preferences]).
    ## This can be done by solving a linear sum assignment problem: the
    ## sought linear order is uniquely characterized by its ranks (r_i),
    ## and the target function is
    ##   \sum_b w_b \sum_i |r_i(b) - r_i| = \sum_{i,k} x_{ik} c_{ik}
    ## where
    ##   c_{ik} = \sum_b w_b | r_i(b) - k |
    ## and x_{ik} is one iff r_i is k.
    ## Clearly, this can be generalized to arbitrary score-based
    ## dissimilarities based on a score function which gives the same
    ## range of values for arbitrary linear orders.

    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")

    nos <- .n_from_control_list(control)

    n <- relation_size(relations)[1L]
    C <- matrix(0, n, n)
    ## Note that Cook and Seiford use Kendall-style "ranks" which are
    ## sorted in decreasing preference, whereas our default scores work
    ## in the opposite direction.
    incidences <- lapply(relations, relation_incidence)
    scores <- sapply(incidences, .incidence_scores_ranks)
    for(k in seq_len(n))
        C[, k] <- rowSums(sweep(abs(scores - k), 2L, weights, "*"))
    .compare <- function(u) outer(u, u, ">=")
    I <- if(nos > 1L)
        lapply(.find_up_to_n_LSAP_solutions(C, nos), .compare)
    else
        .compare(clue::solve_LSAP(C))
    objval <- .relation_consensus_CS_objval(I, incidences, weights)
    meta <- c(.relation_meta_db[["L"]], list(objval = objval))

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

## Consensus methods for central relations using symdiff, Manhattan or
## Euclidean distance.  Note that
##
## * Symdiff dissimilarity only applies to crisp relation, and agrees
##   with Manhattan dissimilarity for these.
##
## * The restricted symdiff fitters are thus only for crisp ensembles.
##
## * We have restricted Manhattan fitters only for crisp ensembles (the
##   symdiff case).  Fitters for fuzzy ensembles could be added by the
##   usual means of turning an l_1 problem with linear/integer
##   constraints into a MILP (see e.g. CLUE).
##
## * The restricted symdiff fitters can also be used for determining
##   restricted *Euclidean* consensus relations for arbitrary (fuzzy)
##   ensembles.  We accomodate for this by calling the internal work
##   horses with a parameter indicating the (symdiff or Euclidean)
##   "context".
##
## * The restricted Euclidean fitters always give crisp relations.
##   Adding fitters for restricted fuzzy consensus would be possible via
##   SUMT approaches (or maybe these could be formulated as mixed
##   integer quadratic programs, but would this help?).

### ** .relation_consensus_Condorcet

## Symdiff/Manhattan crisp consensus relation for an ensemble of crisp
## relations.

.relation_consensus_Condorcet <-
function(relations, weights, control)
{
    if(!.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    nos <- .n_from_control_list(control)

    incidences <- lapply(relations, relation_incidence)
    M <- .make_fit_relation_symdiff_M(incidences, weights)
    I <- M >= 0                         # We do not break ties (>=).
    objval <- .relation_consensus_symdiff_objval(I, incidences, weights)
    meta <- list(is_endorelation = TRUE,
                 is_complete = TRUE,
                 objval = objval)
    if(nos == 1L) {
        meta <- c(meta,
                  list(is_reflexive = all(diag(M) >= 0),
                       is_antisymmetric = all(M != 0)))
        ## According to the way the default solution is defined.
        ## We do not know about transitivity without computing ...
    } else {
        I <- list(I)
        ind <- which(M == 0, arr.ind = TRUE)
        ## Recursively generate up to nos solutions by using 0 or 1 for
        ## the zero entries of M.
        splitter <- function(x, i, j) {
            y <- x
            ## By default we use 1 for the zero entries of M.
            y[i, j] <- 0
            list(x, y)
        }
        k <- 1L
        nr <- nrow(ind)
        while((length(I) < nos) && (k <= nr)) {
            I <- do.call("c",
                         lapply(I, splitter, ind[k, 1L], ind[k, 2L]))
            k <- k + 1L
        }
        if(length(I) > nos)
            I <- I[seq_len(nos)]
    }

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

### ** .relation_consensus_symdiff_A

.relation_consensus_symdiff_A <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "A", weights, control)

### ** .relation_consensus_symdiff_C

.relation_consensus_symdiff_C <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "C", weights, control)

### ** .relation_consensus_symdiff_E

.relation_consensus_symdiff_E <-
function(relations, weights, control)
{
    k <- control$k
    if(!is.null(k))
        .relation_consensus_symdiff_E_k(relations, weights, k, control)
    else
        .relation_consensus_symdiff(relations, "E", weights, control)
}

### ** .relation_consensus_symdiff_L

.relation_consensus_symdiff_L <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "L", weights, control)

### ** .relation_consensus_symdiff_M

.relation_consensus_symdiff_M <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "M", weights, control)

### ** .relation_consensus_symdiff_O

.relation_consensus_symdiff_O <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "O", weights, control)

### ** .relation_consensus_symdiff_W

.relation_consensus_symdiff_W <-
function(relations, weights, control)
{
    k <- control$k
    if(!is.null(k))
        .relation_consensus_symdiff_W_k(relations, weights, k, control)
    else
        .relation_consensus_symdiff(relations, "W", weights, control)
}

### ** .relation_consensus_symdiff_S

.relation_consensus_symdiff_S <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "S", weights, control)

### ** .relation_consensus_symdiff_T

.relation_consensus_symdiff_T <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "T", weights, control)

### ** .relation_consensus_symdiff_preorder

.relation_consensus_symdiff_preorder <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "preorder", weights, control)

### ** .relation_consensus_symdiff_transitive

.relation_consensus_symdiff_transitive <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "transitive", weights, control)

### ** .relation_consensus_manhattan

## Manhattan valued consensus relation for an ensemble of valued
## relations.

.relation_consensus_manhattan <-
function(relations, weights, control)
{
    incidences <- lapply(relations, relation_incidence)
    weights <- rep(weights, length.out = length(incidences))
    ## Incidences of the consensus relation are the weighted medians.
    I <- array(apply(do.call("cbind", lapply(incidences, c)),
                     1L, clue:::weighted_median, weights),
               dim = .size(relations))
    meta <-
        list(objval =
             .relation_consensus_manhattan_objval(I, incidences, weights))
    .make_relation_from_domain_and_incidence(.domain(relations), I, meta)
}

### ** .relation_consensus_euclidean

## Euclidean valued consensus relation for an ensemble of valued
## relations.

.relation_consensus_euclidean <-
function(relations, weights, control)
{
    weights <- rep(weights, length.out = length(relations))
    weights <- weights / sum(weights)
    incidences <- lapply(relations, relation_incidence)
    ## Incidences of the consensus relation are the weighted means.
    I <- .weighted_sum_of_arrays(incidences, weights)
    meta <-
        list(objval =
             .relation_consensus_euclidean_objval(I, incidences, weights))
    .make_relation_from_domain_and_incidence(.domain(relations), I, meta)
}

### ** .relation_consensus_euclidean_A

.relation_consensus_euclidean_A <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "A", weights, control, TRUE)

### ** .relation_consensus_euclidean_C

.relation_consensus_euclidean_C <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "C", weights, control, TRUE)

### ** .relation_consensus_euclidean_E

.relation_consensus_euclidean_E <-
function(relations, weights, control)
{
    k <- control$k
    if(!is.null(k))
        .relation_consensus_symdiff_E_k(relations, weights, k, control, TRUE)
    else
        .relation_consensus_symdiff(relations, "E", weights, control, TRUE)
}

### ** .relation_consensus_euclidean_L

.relation_consensus_euclidean_L <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "L", weights, control, TRUE)

### ** .relation_consensus_euclidean_M

.relation_consensus_euclidean_M <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "M", weights, control, TRUE)

### ** .relation_consensus_euclidean_O

.relation_consensus_euclidean_O <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "O", weights, control, TRUE)

### ** .relation_consensus_euclidean_W

.relation_consensus_euclidean_W <-
function(relations, weights, control)
{
    k <- control$k
    if(!is.null(k))
        .relation_consensus_symdiff_W_k(relations, weights, k, control, TRUE)
    else
        .relation_consensus_symdiff(relations, "W", weights, control, TRUE)
}

### ** .relation_consensus_euclidean_S

.relation_consensus_euclidean_S <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "S", weights, control, TRUE)

### ** .relation_consensus_euclidean_T

.relation_consensus_euclidean_T <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "T", weights, control, TRUE)

### ** .relation_consensus_euclidean_preorder

.relation_consensus_euclidean_preorder <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "preorder",
                                weights, control, TRUE)

### ** .relation_consensus_euclidean_transitive

.relation_consensus_euclidean_transitive <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "transitive",
                                weights, control, TRUE)

## Consensus methods for central relations using CKS distance.
## Currently we only have fitters for crisp ensembles and families where
## the CKS consensus problem can be reduced to a binary linear program.

### ** .relation_consensus_CKS_A

.relation_consensus_CKS_A <-
function(relations, weights, control)
    .relation_consensus_CKS(relations, "A", weights, control)

### ** .relation_consensus_CKS_C

.relation_consensus_CKS_C <-
function(relations, weights, control)
    .relation_consensus_CKS(relations, "C", weights, control)

### ** .relation_consensus_CKS_E

## <NOTE>
## Could add support for CKS/E/k ...
## </NOTE>
.relation_consensus_CKS_E <-
function(relations, weights, control)
    .relation_consensus_CKS(relations, "E", weights, control)

## ** .relation_consensus_CKS_G

.relation_consensus_CKS_G <-
function(relations, weights, control)
{
    ## Generalized (Cook-Kress-Seiford) majority.

    incidences <- lapply(relations, relation_incidence)

    M <- .make_fit_relation_symdiff_M(incidences, weights)
    Q <- .make_fit_relation_CKS_Q(incidences, weights)
    w <- sum(weights)

    I <- (M >= 0) | ((Q < 0) & (M + Q >= - w))
    objval <- .relation_consensus_CKS_objval(I, incidences, weights)
    meta <- list(is_endorelation = TRUE, objval = objval)

    nos <- .n_from_control_list(control)

    if(nos > 1L) {
        I <- list(I)

        ## Split diagonal terms when m_{ii} = 0.
        for(i in which(diag(M) == 0)) {
            I <- do.call("c",
                         lapply(I,
                                function(x, i) {
                                    y <- x
                                    y[i, i] <- 0
                                    list(x, y)
                                },
                                i))
            if(length(I) >= nos) break
        }
        len <- length(I)
        if(len > nos)
            I <- I[seq_len(nos)]
        else if(len < nos) {
            ## Splitting non-diagonal terms is somewhat tricky as we
            ## cannot independently split x_{ij}/x_{ji} for split
            ## positions.
            ## Theory shows that the optimum is characterized via
            ##   \max(q_{ij}, p_{ij}, p_{ji}, 0)
            ## corresponding to 0/0, 1/0, 0/1 and 1/1 for x_{ij}/x_{ji}.
            ## Splits are characterized by
            ##   \max(q_{ij}, p_{ij}, p_{ji}) = 0.
            ## <FIXME>
            ## Check whether the above generalized majority rule always
            ## takes indicidence as one when possible.
            ## </CHECK>
            P <- .make_fit_relation_CKS_P(incidences, weights)
            ## Only need to know where P and Q are zero.
            P <- (P == 0)
            Q <- (Q == 0)
            ## And maybe whether i < j.
            U <- row(P) < col(P)
            ## Helper.
            do_split <- function(I, fun, ind) {
                k <- 1L
                nr <- nrow(ind)
                while((length(I) < nos) && (k <= nr)) {
                    I <- do.call("c",
                                 lapply(I, fun, ind[k, 1L], ind[k, 2L]))
                    k <- k + 1L
                }
                I
            }
            ## If q_{ij} = p_{ij} = 0, can take 0/0 1/0 1/1.
            I <- do_split(I,
                          function(x, i, j) {
                              z <- y <- x
                              x[i, j] <- x[j, i] <- 1
                              y[i, j] <- 1; y[j, i] <- 0
                              z[i, j] <- z[j, i] <- 0
                              list(x, y, z)
                          },
                          which(P & Q, arr.ind = TRUE))
            ## If p_{ij} = p_{ji} = 0, can take 1/0 0/1 1/1.
            I <- do_split(I,
                          function(x, i, j) {
                              z <- y <- x
                              x[i, j] <- x[j, i] <- 1
                              y[i, j] <- 1; y[j, i] <- 0
                              z[i, j] <- 0; z[j, i] <- 1
                              list(x, y, z)
                          },
                          which((P & t(P)) & U, arr.ind = TRUE))
            ## If only q_{ij} = 0, can take 0/1 1/1.
            I <- do_split(I,
                          function(x, i, j) {
                              y <- x
                              x[i, j] <- x[j, i] <- 1
                              y[i, j] <- y[j, i] <- 0
                              list(x, y)
                          },
                          which((Q & !P) & !t(P), arr.ind = TRUE))
            ## If only p_{ij} = 0, can take 1/0 1/1.
            I <- do_split(I,
                          function(x, i, j) {
                              y <- x
                              x[i, j] <- x[j, i] <- 1
                              y[i, j] <- 1; y[j, i] <- 0
                              list(x, y)
                          },
                          which((P & Q) & !t(P), arr.ind = TRUE))
            ## Do this only once and not inside do_split.
            if(length(I) > nos)
                I <- I[seq_len(nos)]
        }
    }

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

### ** .relation_consensus_CKS_L

.relation_consensus_CKS_L <-
function(relations, weights, control)
    .relation_consensus_CKS(relations, "L", weights, control)

### ** .relation_consensus_CKS_M

.relation_consensus_CKS_M <-
function(relations, weights, control)
    .relation_consensus_CKS(relations, "M", weights, control)

### ** .relation_consensus_CKS_O

.relation_consensus_CKS_O <-
function(relations, weights, control)
    .relation_consensus_CKS(relations, "O", weights, control)

### ** .relation_consensus_CKS_W

.relation_consensus_CKS_W <-
function(relations, weights, control)
{
    k <- control$k
    if(!is.null(k))
        .relation_consensus_CKS_W_k(relations, weights, k, control)
    else
        .relation_consensus_CKS(relations, "W", weights, control)
}

### ** .relation_consensus_CKS_S

.relation_consensus_CKS_S <-
function(relations, weights, control)
    .relation_consensus_CKS(relations, "S", weights, control)

### ** .relation_consensus_CKS_T

.relation_consensus_CKS_T <-
function(relations, weights, control)
    .relation_consensus_CKS(relations, "T", weights, control)

### ** .relation_consensus_CKS_preorder

.relation_consensus_CKS_preorder <-
function(relations, weights, control)
    .relation_consensus_CKS_via_QP(relations, "preorder",
                                   weights, control)

### ** .relation_consensus_CKS_transitive

.relation_consensus_CKS_transitive <-
function(relations, weights, control)
    .relation_consensus_CKS_via_QP(relations, "transitive",
                                   weights, control)

### * Relation consensus method registration

## Note that things are simpler here as for CLUE, where we have "typed"
## db's (partition or hierarchy type).

## We currently do without explicit registry getters and setters, which
## in the non-typed case could simplify to the following:
##   get_methods_from_db <-
##   function(db)
##       objects(db)
##   get_method_from_db <-
##   function(db, name)
##       db[[name]]
##   put_method_into_db <-
##   function(db, name, value)
##       db[[name]] <- value
## However, we provide a getter which allows for partial matching.

relation_consensus_methods_db <- new.env()
get_relation_consensus_method <-
function(name)
{
    keys <- objects(relation_consensus_methods_db)
    ind <- pmatch(name, keys)
    if(is.na(ind))
        stop(gettextf("Invalid consensus method '%s'.", name),
             domain = NA)
    relation_consensus_methods_db[[keys[ind]]]
}
set_relation_consensus_method <-
function(name, definition, ...)
{
    ## Note that consensus methods are not necessarily optimization
    ## based (and hence do not necessarily have associated dissimilarity
    ## and exponent).
    value <- c(list(definition = definition), list(...))
    class(value) <- "relation_consensus_method"
    relation_consensus_methods_db[[name]] <- value
}

set_relation_consensus_method("Borda",
                              .relation_consensus_Borda)
set_relation_consensus_method("Copeland",
                              .relation_consensus_Copeland)
set_relation_consensus_method("majority",
                              .relation_consensus_majority)
## Note that constructive methods do not necessarily give central
## relations.
set_relation_consensus_method("CKS/A",
                              .relation_consensus_CKS_A,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/C",
                              .relation_consensus_CKS_C,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/E",
                              .relation_consensus_CKS_E,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/G",
                              .relation_consensus_CKS_G,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/L",
                              .relation_consensus_CKS_L,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/M",
                              .relation_consensus_CKS_M,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/O",
                              .relation_consensus_CKS_O,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/S",
                              .relation_consensus_CKS_S,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/T",
                              .relation_consensus_CKS_T,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/W",
                              .relation_consensus_CKS_W,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/preorder",
                              .relation_consensus_CKS_preorder,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/transitive",
                              .relation_consensus_CKS_transitive,
                              dissimilarity = "CKS",
                              exponent = 1)

set_relation_consensus_method("CS",
                              .relation_consensus_CS,
                              dissimilarity = "CS",
                              exponent = 1)
set_relation_consensus_method("Condorcet",
                              .relation_consensus_Condorcet,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/A",
                              .relation_consensus_symdiff_A,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/C",
                              .relation_consensus_symdiff_C,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/E",
                              .relation_consensus_symdiff_E,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/G",
                              .relation_consensus_Condorcet,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/L",
                              .relation_consensus_symdiff_L,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/M",
                              .relation_consensus_symdiff_M,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/O",
                              .relation_consensus_symdiff_O,
                              dissimilarity = "symdiff",
                              exponent = 1)
## <NOTE>
## Keep this for back-compatibility.
set_relation_consensus_method("SD/P",
                              .relation_consensus_symdiff_W,
                              dissimilarity = "symdiff",
                              exponent = 1)
## </NOTE>
set_relation_consensus_method("SD/S",
                              .relation_consensus_symdiff_S,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/T",
                              .relation_consensus_symdiff_T,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/W",
                              .relation_consensus_symdiff_W,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/preorder",
                              .relation_consensus_symdiff_preorder,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/transitive",
                              .relation_consensus_symdiff_transitive,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("manhattan",
                              .relation_consensus_manhattan,
                              dissimilarity = "manhattan",
                              exponent = 1)
set_relation_consensus_method("euclidean",
                              .relation_consensus_euclidean,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/A",
                              .relation_consensus_euclidean_A,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/C",
                              .relation_consensus_euclidean_C,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/E",
                              .relation_consensus_euclidean_E,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/G",
                              .relation_consensus_Condorcet,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/L",
                              .relation_consensus_euclidean_L,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/M",
                              .relation_consensus_euclidean_M,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/O",
                              .relation_consensus_euclidean_O,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/S",
                              .relation_consensus_euclidean_S,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/T",
                              .relation_consensus_euclidean_T,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/W",
                              .relation_consensus_euclidean_W,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/preorder",
                              .relation_consensus_euclidean_preorder,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/transitive",
                              .relation_consensus_euclidean_transitive,
                              dissimilarity = "euclidean",
                              exponent = 2)

### * Relation consensus workers

### ** .relation_consensus_symdiff

.relation_consensus_symdiff <-
function(relations, family, weights, control, euclidean = FALSE)
{
    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")
    if(!euclidean && !.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    incidences <- lapply(relations, relation_incidence)
    M <- .make_fit_relation_symdiff_M(incidences, weights)
    I <- fit_relation_symdiff(M, family, control)
    objval <- if(euclidean)
        .relation_consensus_euclidean_objval(I, incidences, weights)
    else
        .relation_consensus_symdiff_objval(I, incidences, weights)
    meta <- c(.relation_meta_db[[family]], list(objval = objval))

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

### ** .relation_consensus_symdiff_E_k

.relation_consensus_symdiff_E_k <-
function(relations, weights, k, control, euclidean = FALSE)
{
    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")
    if(!euclidean && !.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    incidences <- lapply(relations, relation_incidence)
    M <- .make_fit_relation_symdiff_M(incidences, weights)
    I <- fit_relation_symdiff_E_k(M, k, control)
    objval <- if(euclidean)
        .relation_consensus_euclidean_objval(I, incidences, weights)
    else
        .relation_consensus_symdiff_objval(I, incidences, weights)
    meta <- c(.relation_meta_db[["E"]], list(objval = objval))

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

### ** .relation_consensus_symdiff_W_k

.relation_consensus_symdiff_W_k <-
function(relations, weights, k, control, euclidean = FALSE)
{
    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")
    if(!euclidean && !.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    incidences <- lapply(relations, relation_incidence)
    M <- .make_fit_relation_symdiff_M(incidences, weights)
    I <- fit_relation_symdiff_W_k(M, k, control)
    objval <- if(euclidean)
        .relation_consensus_euclidean_objval(I, incidences, weights)
    else
        .relation_consensus_symdiff_objval(I, incidences, weights)
    meta <- c(.relation_meta_db[["W"]], list(objval = objval))

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

## ** .relation_consensus_CKS

.relation_consensus_CKS <-
function(relations, family, weights, control)
{
    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")
    if(!.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    ## Solve via QP formulation if requested explicitly.
    ## Mostly useful for testing purposes.
    if(identical(control$simplify, FALSE))
       return(.relation_consensus_CKS_via_QP(relations, family,
                                             weights, control))

    incidences <- lapply(relations, relation_incidence)
    if(family %in% c("C", "L", "M", "T", "W")) {
        ## For the complete families, we can fit CKS via transforming
        ## incidences and fitting symdiff.
        incidences <-
            lapply(incidences,
                   function(I) {
                       ind <- row(I) != col(I)
                       I[ind] <- pmax(I, 1 - t(I))[ind]
                       I
                   })
        M <- .make_fit_relation_symdiff_M(incidences, weights)
    } else if(family %in% c("O", "A")) {
        ## If the consensus relation is antisymmetric, the CKS consensus
        ## problem reduces to a linear program with coefficients
        ## $p_{ij} - q_{ij}$ and $- q_{ii}$ for the off-diagonal and
        ## diagonal terms, respectively.
        ## For O, the diagonal is "fixed" by the family.
        ## For A, it is handled by the majority rule.
        P <- .make_fit_relation_CKS_P(incidences, weights)
        Q <- .make_fit_relation_CKS_Q(incidences, weights)
        M <- P - Q
        ## Note that $p_{ii} = \sum_b w_b (1 - 2 p_{ii}(b)) = \sum_b w_b.
        diag(M) <- diag(M) - sum(weights)
        ## Alternatively, do diag(M) <- - diag(Q).
    } else if(family %in% c("E", "S")) {
        ## If the consensus relation is symmetric, the CKS consensus
        ## problem reduces to a linear program with coefficients
        ## $- q_{ij} / 2$ and $- q_{ii}$ for the off-diagonal and
        ## diagonal terms, respectively.
        ## For E, the diagonal is "fixed" by the family.
        ## For S, it is handled by the majority diagonal rule
        ## Note that it is not necessary to halve the off-diagonal terms
        ## as the diagonal and off-diagonal terms are handled separately
        ## by symdiff fitting.
        M <- - .make_fit_relation_CKS_Q(incidences, weights)
    } else {
        stop("Not implemented.")
    }
    I <- fit_relation_symdiff(M, family, control)
    objval <- .relation_consensus_CKS_objval(I, incidences, weights)
    meta <- c(.relation_meta_db[[family]], list(objval = objval))

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

## ** .relation_consensus_CKS_W_k

.relation_consensus_CKS_W_k <-
function(relations, weights, k, control)
{
    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")
    if(!.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    incidences <-
        lapply(relations,
               function(R) {
                   I <- relation_incidence(R)
                   ind <- row(I) != col(I)
                   I[ind] <- pmax(I, 1 - t(I))[ind]
                   I
               })
    M <- .make_fit_relation_symdiff_M(incidences, weights)
    I <- fit_relation_symdiff_W_k(M, k, control)
    objval <- .relation_consensus_CKS_objval(I, incidences, weights)
    meta <- c(.relation_meta_db[["W"]], list(objval = objval))

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

## ** .relation_consensus_CKS_via_QP

.relation_consensus_CKS_via_QP <-
function(relations, family, weights, control)
{
    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")
    if(!.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    incidences <- lapply(relations, relation_incidence)
    P <- .make_fit_relation_CKS_P(incidences, weights)
    Q <- .make_fit_relation_CKS_Q(incidences, weights)
    I <- fit_relation_CKS_via_QP(P, Q, family, control)
    objval <- .relation_consensus_CKS_objval(I, incidences, weights)
    meta <- c(.relation_meta_db[[family]], list(objval = objval))

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

### * Utilities

### ** .make_fit_relation_symdiff_M

.make_fit_relation_symdiff_M <-
function(incidences, weights)
{
    ## Compute the array
    ##   \sum_b w_b (2 incidence(b) - 1)
    ## used in fit_relation_symdiff() and also for the Condorcet
    ## consensus solution.

    w <- rep(weights, length.out = length(incidences))
    if(is.relation_ensemble(incidences))
        incidences <- lapply(incidences, relation_incidence)
    M <- .weighted_sum_of_arrays(incidences, w)
    2 * M - sum(w)
}

### ** .make_fit_relation_CKS_P

.make_fit_relation_CKS_P <-
function(incidences, weights)
    .make_fit_relation_symdiff_M(lapply(incidences,
                                        function(I) pmin(I, 1 - t(I))),
                                 weights)

### ** .make_fit_relation_CKS_Q

.make_fit_relation_CKS_Q <-
function(incidences, weights)
    .make_fit_relation_symdiff_M(lapply(incidences,
                                        function(I) 1 - pmax(I, t(I))),
                                 weights)


### ** .make_consensus_from_incidences

.make_consensus_from_incidences <-
function(D, I, meta = NULL)
{
    if(is.list(I)) {
        ## In case *all* consensus solutions were sought ...
        relations <-
            lapply(I,
                   function(e)
                   .make_relation_from_domain_and_incidence(D, e, meta))
        relation_ensemble(list = relations)
    }
    else
        .make_relation_from_domain_and_incidence(D, I, meta)
}

## Utilities for computing the values of the objective functions for the
## optimization-based consensus methods.

### ** .relation_consensus_symdiff_objval

.relation_consensus_symdiff_objval <-
function(I, incidences, weights)
{
    ## Be nice.
    if(is.relation_ensemble(incidences))
        incidences <- lapply(incidences, relation_incidence)
    if(is.list(I)) I <- I[[1L]]

    sum(weights *
        sapply(incidences, .incidence_dissimilarity_symdiff, I))
}

### ** .relation_consensus_manhattan_objval

.relation_consensus_manhattan_objval <-
function(I, incidences, weights)
{
    ## Be nice.
    if(is.relation_ensemble(incidences))
        incidences <- lapply(incidences, relation_incidence)
    if(is.list(I)) I <- I[[1L]]

    sum(weights *
        sapply(incidences, .incidence_dissimilarity_manhattan, I))
}

### ** .relation_consensus_euclidean_objval

.relation_consensus_euclidean_objval <-
function(I, incidences, weights)
{
    ## Be nice.
    if(is.relation_ensemble(incidences))
        incidences <- lapply(incidences, relation_incidence)
    if(is.list(I)) I <- I[[1L]]

    sum(weights *
        sapply(incidences, .incidence_dissimilarity_euclidean, I) ^ 2)
}

### ** .relation_consensus_CS_objval

.relation_consensus_CS_objval <-
function(I, incidences, weights)
{
    ## Be nice.
    if(is.relation_ensemble(incidences))
        incidences <- lapply(incidences, relation_incidence)
    if(is.list(I)) I <- I[[1L]]

    sum(weights * sapply(incidences, .incidence_dissimilarity_CS, I))
}

### ** .relation_consensus_CKS_objval

.relation_consensus_CKS_objval <-
function(I, incidences, weights)
{
    ## Be nice.
    if(is.relation_ensemble(incidences))
        incidences <- lapply(incidences, relation_incidence)
    if(is.list(I)) I <- I[[1L]]

    sum(weights * sapply(incidences, .incidence_dissimilarity_CKS, I))
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
