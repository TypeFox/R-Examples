## Relation fitters.

### * fit_relation_symdiff

## General purpose fitter for families
.SD_families <- c("E", "L", "O", "W", "T", "C", "A", "S", "M",
                  "preorder", "transitive")
## where:
##   E ............ Equivalence relations                 REF SYM TRA
##   L ............ Linear orders                     TOT REF ASY TRA
##   O ............ Partial orders                        REF ASY TRA
##   W ............ Weak orders (complete preorders)  TOT REF     TRA
##   T ............ Tournament                        TOT IRR ASY
##   C ............ Complete relations                TOT
##   A ............ Antisymmetric relations                   ASY
##   S ............ Symmetric relations                       SYM
##   M ............ Matches                           TOT REF
##   preorder ..... Preorders                             REF     TRA
##   transitive ... Transitive relations                          TRA
## and
##   TOT ... total/complete
##   REF ... reflexive
##   ASY ... antisymmetric,
##   IRR ... irreflexive
##   TRA ... transitive
##   SYM ... symmetric
##
## Families which are not REF or IRR by definition will be reflexive iff
## the (weighted) majority of all input relations is.

## Keep this is sync with .relation_meta_db (currently in utilities.R).

## <NOTE>
## Ideally we would move to a database for families which specify their
## parametrizations (currently, upper triangular or off-diagonal) and
## constraints etc., in particular making it more convenient and
## reliable to add support for new families.
## A single classification (e.g., upper_tri vs. offdiag) does not do the
## job.
## </NOTE>

## Number of non-redundant object pairs.
.n_of_pairs <-
function(n, family)
{
    family <- match.arg(family, .SD_families)
    N <- n * (n - 1L)                   # Number of distinct pairs.
    switch(EXPR = family,
           E =, L =, S =, T = N / 2,    # upper_tri parametrization
           A =, C =, M =, O =, W =,
           preorder =, transitive = N   # offdiag parametrization
           )
}

## Note that for families which are symmetric (x_{ij} = x_{ji}, E and S)
## or complete and antisymmetric (x_{ij} = 1 - x_{ji}, L and T) we use
## an upper_tri parametrization.  For the others, we use offdiag:
.SD_families_using_offdiag_parametrization <-
    .SD_families[!.SD_families %in% c("E", "L", "S", "T")]
## Some of these have TOT or ASY constraints (but not both):
.SD_families_using_tot_or_asy_constraints <-
    c("A", "C", "M", "O", "W")
## Note that for the families using upper_tri parametrization, TOT or
## ASY constraints either give trivial (together with SYM: E and S)
## results or are redundant.

## Number of transitivity constraints for the non-redundant pairs.
## (None for tournaments.)
.n_of_transitivity_constraints <-
function(n, family)
{
    family <- match.arg(family, .SD_families)
    N <- n * (n - 1L) * (n - 2L)        # Number of distinct triples.
    switch(EXPR = family,
           E = N / 2,
           L = N / 3,
           O =, W =, preorder =, transitive = N,
           A =, C =, M =, S =, T = 0L   # Families w/out TRA.
           )
}

## Make function giving the position of incidence (i, j) in the vector
## of non-redundant incidences used for symdiff fitting (upper.tri() for
## E/L/S/T and .offdiag() for the others, respectively).
.make_pos <-
function(n, family)
{
    family <- match.arg(family, .SD_families)
    if(family %in% .SD_families_using_offdiag_parametrization)
        function(i, j) {
            ## Position of x_{ij} in x[.offdiag(x)]
            (n - 1L) * (j - 1L) + i - (i >= j)
        }
    else
        function(i, j) {
            ## Position of x_{ij} in x[upper.tri(x)].
            i + (j - 1L) * (j - 2L) / 2
        }
}

fit_relation_symdiff <-
function(x, family = .SD_families, control = list())
{
    sparse <- !identical(control$sparse, FALSE)

    ## Number of objects:
    n <- nrow(x)
    objective_in <- if (family %in% c("L", "T"))
        (x - t(x))[upper.tri(x)]
    else if (family %in% c("E", "S"))
        (x + t(x))[upper.tri(x)]
    else
        x[.offdiag(x)]

    ## Handle constraints implied by the family to be fitted.
    ## Need all variables in { 0 , 1 }, i.e., >= 0 and <= 1 and binary,
    ## plus maybe totality or antisymmetry, plus maybe transitivity.
    NP <- .n_of_pairs(n, family)
    eye <-
        if(!sparse) diag(1, NP) else simple_triplet_diag_matrix(1, NP)
    constr_mat <-
        rbind(eye,
              eye,
              if(family %in% .SD_families_using_tot_or_asy_constraints)
                  .make_tot_or_asy_constraint_mat(n, sparse),
              .make_transitivity_constraint_mat(n, family, sparse))
    constr_dir <-
        c(rep.int(">=", NP),
          rep.int("<=", NP),
          if(family %in% .SD_families_using_tot_or_asy_constraints)
              .make_tot_or_asy_constraint_dir(n, family),
          .make_transitivity_constraint_dir(n, family))
    constr_rhs <-
        c(rep.int(0, NP),
          rep.int(1, NP),
          if(family %in% .SD_families_using_tot_or_asy_constraints)
              .make_tot_or_asy_constraint_rhs(n),
          .make_transitivity_constraint_rhs(n, family))

    ## Handle additional constraints.
    acmaker <-
        .make_additional_constraint_maker_using_incidences(n, family, sparse)
    if(!is.null(A <- control$constraints)) {
        A <- .canonicalize_additional_constraints(A, family)
        add <- acmaker(A)
        constr_mat <- rbind(constr_mat, add$mat)
        constr_dir <- c(constr_dir, add$dir)
        constr_rhs <- c(constr_rhs, add$rhs)
    }

    labels <- dimnames(x)

    nos <- .n_from_control_list(control)
    one <- nos == 1L

    ## Compute diagonal:
    ## For families which are reflexive or irreflexive, set all ones or
    ## zeroes, respectively.  For all other families, a diagonal entry
    ## will be set iff it is in the (non-strict weighted) majority of
    ## the input relations.
    diagonal <- if(family %in% c("E", "L", "M", "O", "W", "preorder"))
        rep.int(1, n)
    else if(family == "T")
        rep.int(0, n)
    else if(one)
        diag(x) >= 0
    else {
        ## If more than one solution is sought, mark diagonal entries
        ## with no strict majority as NA to allow for later expansion.
        ifelse(diag(x) != 0, diag(x) >= 0, NA)
    }

    types <- rep.int("B", NP)
    milp <- MILP(objective_in,
                 list(constr_mat,
                      constr_dir,
                      constr_rhs),
                 types = types,
                 maximum = TRUE)
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")
    
    if(!one) {
        ## Do this in a separate function which also handles diagonal
        ## expansion.
        return(.find_up_to_n_relation_symdiff_optima(milp,
                                                     nos,
                                                     control$solver,
                                                     control$control,
                                                     family, labels,
                                                     diagonal,
                                                     verbose))
    }
    out <- solve_MILP(milp, control$solver,
                      c(list(verbose = verbose), control$control))
    .stop_if_lp_status_is_nonzero(out$status, family)
    ## Turn the solution back into a full incidence matrix.
    fit <- .make_incidence_from_symdiff_MILP_solution(out$solution,
                                                      family, labels,
                                                      diagonal)
    ## For the time being, tack some of the MILP results on so that we
    ## can look at them (but not everything due to size ...)
    ## <FIXME>
    ## Remove this eventually ...
    attr(fit, ".lp") <- out[c("solution", "objval")]
    ## </FIXME>
    fit
}

### * .stop_if_lp_status_is_nonzero

.stop_if_lp_status_is_nonzero <-
function(status, family)
{
    if(status != 0) {
        ## This should really only be possible in case additional
        ## constraints were given.
        stop(gettextf("Given constraints are incompatible with family '%s'.",
                      family),
             domain = NA)
    }
}

### * .find_up_to_n_relation_symdiff_optima

.find_up_to_n_relation_symdiff_optima <-
function(x, n, solver, control, family, labels, diagonal, verbose)
{
    ## Note that this only gets used if more than one solution is
    ## sought.
    y <- solve_MILP(x, solver,
                    c(list(n = n, verbose = verbose), control))
    ## Check status:
    if(length(y) == 1L)
        .stop_if_lp_status_is_nonzero(y[[1L]]$status, family)
    ## Turn back into solutions, expanding diagonal entries if needed.
    .make_incidence <- function(e, d)
        .make_incidence_from_symdiff_MILP_solution(e$solution,
                                                   family, labels, d)
    y <- if(!any(ind <- which(is.na(diagonal))))
        lapply(y, .make_incidence, diagonal)
    else {
        diagonals <- list(diagonal)
        for(i in ind) {
            diagonals <- do.call(c,
                                 lapply(diagonals,
                                        function(x) {
                                            u <- v <- x
                                            u[i] <- 0
                                            v[i] <- 1
                                            list(u, v)
                                        }))
        }
        do.call(c,
                lapply(diagonals,
                       function(diagonal)
                       lapply(y, .make_incidence, diagonal)))
    }
    if(length(y) > n)
        y <- y[seq_len(n)]
    y
}    

### * Constraint generators: transitivity.

### ** .make_transitivity_constraint_mat

.make_transitivity_constraint_mat <-
function(n, family, sparse = FALSE)
{
    family <- match.arg(family, .SD_families)

    NP <- .n_of_pairs(n, family)
    if ((n <= 2L) || (family %in% c("A", "C", "M", "S", "T"))) {
        if(!sparse)
            return(matrix(0, 0, NP))
        else
            return(simple_triplet_zero_matrix(0, NP))
    }

    NC <- .n_of_transitivity_constraints(n, family)
    pos <- .make_pos(n, family)

    if(family %in% c("E", "L")) {
        ## Create a matrix with all combinations of triples i < j < k in
        ## the rows.
        ind <- seq_len(n)
        z <- as.matrix(expand.grid(ind, ind, ind))[, c(3L, 2L, 1L)]
        z <- z[(z[, 1L] < z[, 2L]) & (z[, 2L] < z[, 3L]), , drop = FALSE]

        p_ij <- pos(z[, 1L], z[, 2L])
        p_ik <- pos(z[, 1L], z[, 3L])
        p_jk <- pos(z[, 2L], z[, 3L])

        ind <- seq_len(NC)

        if(!sparse) {
            out <- matrix(0, NC, NP)
            if(family == "E") {
                ## For equivalence relations, we have 3 transitivity
                ## constraints for each such triple.
                out[cbind(ind, c(p_ij, p_ij, p_ik))] <- 1
                out[cbind(ind, c(p_jk, p_ik, p_jk))] <- 1
                out[cbind(ind, c(p_ik, p_jk, p_ij))] <- -1
            }
            else if(family == "L") {
                ## For linear orders, we only get two constraints.
                NT <- NC / 2
                out[cbind(ind, c(p_ij, p_ij))] <-
                    rep(c(1, -1), each = NT)
                out[cbind(ind, c(p_jk, p_jk))] <-
                    rep(c(1, -1), each = NT)
                out[cbind(ind, c(p_ik, p_ik))] <-
                    rep(c(-1, 1), each = NT)
            }
        } else {
            out <- if(family == "E")
                simple_triplet_matrix(rep.int(ind, 3L),
                                      c(p_ij, p_ij, p_ik,
                                        p_jk, p_ik, p_jk,
                                        p_ik, p_jk, p_ij),
                                      c(rep.int(1, 2L * NC),
                                        rep.int(-1, NC)),
                                      NC, NP)
            else if(family == "L") {
                NT <- NC / 2
                simple_triplet_matrix(rep.int(ind, 3L),
                                      c(p_ij, p_ij,
                                        p_jk, p_jk,
                                        p_ik, p_ik),
                                      c(rep(c(1, -1), each = NT),
                                        rep(c(1, -1), each = NT),
                                        rep(c(-1, 1), each = NT)),
                                      NC, NP)
            }
        }
    }
    else {
        ## Create a matrix with all combinations of distinct triples i,
        ## j, k in the rows.
        ind <- seq_len(n)
        z <- as.matrix(expand.grid(ind, ind, ind))[, c(3L, 2L, 1L)]
        z <- z[(z[, 1L] != z[, 2L])
               & (z[, 2L] != z[, 3L])
               & (z[, 3L] != z[, 1L]), ]

        ind <- seq_len(NC)
        if(!sparse) {
            out <- matrix(0, NC, NP)
            out[cbind(ind, pos(z[, 1L], z[, 2L]))] <- 1
            out[cbind(ind, pos(z[, 2L], z[, 3L]))] <- 1
            out[cbind(ind, pos(z[, 1L], z[, 3L]))] <- -1
        } else {
            out <- simple_triplet_matrix(rep.int(ind, 3L),
                                         c(pos(z[, 1L], z[, 2L]),
                                           pos(z[, 2L], z[, 3L]),
                                           pos(z[, 1L], z[, 3L])),
                                         c(rep.int(1, 2L * NC),
                                           rep.int(-1, NC)),
                                         NC, NP)
        }
    }

    out
}

### ** .make_transitivity_constraint_dir

.make_transitivity_constraint_dir <-
function(n, family)
{
    if(n <= 2L) return(character())
    family <- match.arg(family, .SD_families)
    rep.int("<=", .n_of_transitivity_constraints(n, family))
}

### ** .make_transitivity_constraint_rhs

.make_transitivity_constraint_rhs <-
function(n, family)
{
    if(n <= 2L) return(double())

    family <- match.arg(family, .SD_families)

    NC <- .n_of_transitivity_constraints(n, family)

    if(family == "L")
        rep(c(1, 0), each = NC / 2)
    else
        rep.int(1, NC)
}

### * Constraint generators: completeness or antisymmetry.

## This translates into
##
##    x_{ij} + x_{ji} >= 1      [completeness: C M W]
##    x_{ij} + x_{ji} <= 1      [antisymmetry: A O]
##
## For tournaments, we have both, and hence x_{ij} + x_{ji} = 1 as for
## linear orders, such that we use the non-redundant upper diagonal
## pairs representation, and get no additional explicit constraints.
##
## As these constraints only differ by direction, we handle them
## together.

### ** .make_tot_or_asy_constraint_mat

.make_tot_or_asy_constraint_mat <-
function(n, sparse = FALSE)
{
    if(n <= 1L) {
        if(!sparse)
            return(matrix(0, 0L, 0L))   # :-)
        else
            return(simple_triplet_zero_matrix(0L, 0L))
    }

    ## Position of x_{ij} in x[.offdiag(x)]
    pos <- .make_pos(n, "W")

    ## Number of non-redundant pairs.
    NP <- n * (n - 1L)
    ## Number of constraints.
    NC <- NP / 2

    ind <- seq_len(n)
    z <- as.matrix(expand.grid(ind, ind))[, c(2L, 1L)]
    z <- z[z[, 1L] < z[, 2L], , drop = FALSE]

    ind <- seq_len(NC)

    if(!sparse) {
        out <- matrix(0, NC, NP)
        out[cbind(ind, pos(z[, 1L], z[, 2L]))] <- 1
        out[cbind(ind, pos(z[, 2L], z[, 1L]))] <- 1
    }
    else {
        out <- simple_triplet_matrix(c(ind, ind),
                                     c(pos(z[, 1L], z[, 2L]),
                                       pos(z[, 2L], z[, 1L])),
                                     rep.int(1, 2L * NC),
                                     NC, NP)
    }

    out
}

### ** .make_tot_or_asy_constraint_dir

.make_tot_or_asy_constraint_dir <-
function(n, family)
{
    NC <- n * (n - 1L) / 2
    rep.int(switch(EXPR = family,
                   A =, O = "<=",       # ASY
                   C =, M =, W = ">="   # TOT
                   ),
            NC)
}

### ** .make_tot_or_asy_constraint_rhs

.make_tot_or_asy_constraint_rhs <-
function(n)
{
    NC <- n * (n - 1L) / 2
    rep.int(1, NC)
}

### * Constraint generators: explicitly given constraints.

### .make_additional_constraint_maker_using_incidences

.make_additional_constraint_maker_using_incidences <-
function(n, family, sparse = FALSE)
    function(A) {
        ## (This assumes A has already been validated and
        ## canonicalized.)
        na <- nrow(A)
        nc <- .n_of_pairs(n, family)
        pos <- .make_pos(n, family)
        if(!sparse) {
            mat <- matrix(0, na, nc)
            mat[cbind(seq_len(na), pos(A[, 1L], A[, 2L]))] <- 1
        } else {
            mat <- simple_triplet_matrix(seq_len(na),
                                         pos(A[, 1L], A[, 2L]),
                                         rep.int(1, na),
                                         na, nc)
        }
        list(mat = mat,
             dir = rep.int("==", na),
             rhs = A[, 3L])
    }

### .make_additional_constraint_maker_using_memberships_E

.make_additional_constraint_maker_using_memberships_E <-
function(pos, n_of_variables, nc, sparse = FALSE)
    function(A) {
        ## (This assumes A has already been validated and
        ## canonicalized.)
        na <- nrow(A)
        ## For equivalences (E):
        ## A constraint of the form (i, j, 1) implies that
        ## objects i and j are in the same class, i.e.:
        ##   m_{ik} - m_{jk} == 0   for all k.
        ## A constraint of the form (i, j, 0) implies that
        ## objects i and j are in different classes, i.e.:
        ##   m_{ik} + m_{jk} <= 1   for all k.
        ind <- seq_len(na)
        if(!sparse) {
            mat <- matrix(0, nc * na, n_of_variables)
            for(k in seq_len(nc)) {
                mat[cbind(ind, pos(A[, 1L], k))] <- 1
                mat[cbind(ind, pos(A[, 2L], k))] <-
                    ifelse(A[, 3L] == 1, -1, 1)
                ind <- ind + na
            }
        } else {
            i <- rep.int(ind, 2L * nc) +
                rep(seq(from = 0, by = na, length.out = nc),
                    each = 2L * na)
            j <- unlist(lapply(seq_len(nc),
                               function(k)
                               c(pos(A[, 1L], k), pos(A[, 2L], k))))
            v <- rep.int(c(rep.int(1, na),
                           ifelse(A[, 3L] == 1, -1, 1)),
                         nc)
            mat <- simple_triplet_matrix(i, j, v,
                                         nc * na, n_of_variables)
        }
        list(mat = mat,
             dir = rep.int(ifelse(A[, 3L] == 1, "==", "<="), nc),
             rhs = rep.int(ifelse(A[, 3L] == 1, 0, 1), nc))
    }

### .make_additional_constraint_maker_using_memberships_W

.make_additional_constraint_maker_using_memberships_W <-
function(pos, n_of_variables, nc, sparse = FALSE)
    function(A) {
        ## (This assumes A has already been validated and
        ## canonicalized.)
        na <- nrow(A)
        ## For weak orders (W):
        ## A constraint of the form (i, j, 1) means that i <= j,
        ## or class(i) <= class(j).
        ## I.e., can only have m_{jk} = 1 if \sum_{l <= k} m_{il} = 1,
        ## giving
        ##   m_{jk} <= \sum_{l <= k} m_{il}   for all k.
        ## A constraint of the form (i, j, 0) means that i > j,
        ## or class(i) > class(j).
        ## I.e., can only have m_{ik} = 1 if \sum_{l < k} m_{jl} = 1,
        ## giving
        ##   m_{ik} <= \sum_{l < k} m_{jl}    for all k,
        ## (with the empty sum for k = 1 taken as zero).
        if(!sparse) {
            mat <- matrix(0, nc * na, n_of_variables)
            ind <- i1 <- which(A[, 3L] == 1)
            if(any(i1)) {
                for(k in seq_len(nc)) {
                    mat[cbind(ind, pos(A[i1, 2L], k))] <- 1
                    for(l in seq_len(k))
                        mat[cbind(ind, pos(A[i1, 1L], l))] <- -1
                    ind <- ind + na
                }
            }
            ind <- i0 <- which(A[, 3L] == 0)
            if(any(i0)) {
                for(k in seq_len(nc)) {
                    mat[cbind(ind, pos(A[i0, 1L], k))] <- 1
                    for(l in seq_len(k - 1L))
                        mat[cbind(ind, pos(A[i0, 2L], l))] <- -1
                    ind <- ind + na
                }
            }
        } else {
            i0 <- i1 <- j0 <- j1 <- integer()
            v0 <- v1 <- double()
            ## This is somewhat messy to compute explicitly in one pass,
            ## so let's simply use a loop for building things up (but
            ## avoid looping over j <- c(j, something) constructs for
            ## performance reasons.
            ind <- which(A[, 3L] == 0)
            if(any(ind)) {
                len <- length(ind)
                i0 <- unlist(lapply(seq_len(nc),
                                    function(k)
                                    rep.int(ind + (k - 1L) * na, k)
                                    ))
                j0 <- unlist(lapply(seq_len(nc),
                                    function(k)
                                    c(pos(A[ind, 1L], k),
                                      unlist(lapply(seq_len(k - 1L),
                                                    function(l)
                                                    pos(A[ind, 2L], l))))
                                    ))
                v0 <- unlist(lapply(seq_len(nc),
                                    function(k)
                                    c(rep.int(1, len),
                                      rep.int(-1, (k - 1L) * len))))
            }
            ind <- which(A[, 3L] == 1)
            if(any(ind)) {
                len <- length(ind)
                i1 <- unlist(lapply(seq_len(nc),
                                    function(k)
                                    rep.int(ind + (k - 1L) * na, k + 1L)
                                    ))
                j1 <- unlist(lapply(seq_len(nc),
                                    function(k)
                                    c(pos(A[ind, 2L], k),
                                      unlist(lapply(seq_len(k),
                                                    function(l)
                                                    pos(A[ind, 1L], l))))
                                    ))
                v1 <- unlist(lapply(seq_len(nc),
                                    function(k)
                                    c(rep.int(1, len),
                                      rep.int(-1, k * len))))
            }
            mat <- simple_triplet_matrix(c(i0, i1),
                                         c(j0, j1),
                                         c(v0, v1),
                                         nc * na, n_of_variables)
        }

        list(mat = mat,
             dir = rep.int("<=", nc * na),
             rhs = double(nc * na))
    }

### ** .make_additional_constraint_maker_using_o_c_pairs

.make_additional_constraint_maker_using_o_c_pairs <-
function(pos, n_of_variables, nc, sparse = FALSE)
    function(A) {
        ## (This assumes A has already been validated and
        ## canonicalized.)
        na <- nrow(A)
        if(!sparse) {
            mat <- matrix(0, na, n_of_variables)
            mat[cbind(seq_len(na), pos(A[, 1L], A[, 2L]))] <- 1
        } else {
            mat <- simple_triplet_matrix(seq_len(na),
                                         pos(A[, 1L], A[, 2L]),
                                         rep.int(1, na),
                                         na, n_of_variables)
        }
        list(mat = mat, dir = rep.int("==", na), rhs = A[, 3L])
    }

### ** .canonicalize_additional_constraints

.canonicalize_additional_constraints <-
function(A, family)
{
    ## Additional constraints should be given as a 3-column matrix with
    ## rows (i, j, x) meaning that the incidence of i and j should be
    ## equal to x.
    if(!is.matrix(A) || (ncol(A) != 3L))
        stop("Invalid additional incidence constraints.")
    ## For all families, we only use off-diagonal terms, so drop
    ## diagonal ones (which should all be one, of course).
    A <- A[A[, 1L] != A[, 2L], , drop = FALSE]
    if(!nrow(A)) return(matrix(0, 0L, 3L))
    ## For E/L/S/T, we use only pairs i < j, so swap and possibly flip
    ## (L/T) if needed.
    if(family %in% c("E", "L", "S", "T")) {
        ind <- A[, 1L] > A[, 2L]
        if(any(ind))
            A[ind, ] <- cbind(A[ind, c(2L, 1L), drop = FALSE],
                              if(family %in% c("E", "S")) A[ind, 3L]
                              else 1 - A[ind, 3L])
    }
    ## Now validate.
    pos <- max(A) * (A[, 2L] - 1) + A[, 1L]
    ind <- duplicated(pos)
    if(any(ind)) {
        ## Check if the duplicated entries all have the same
        ## incidences.
        lens <- tapply(A[, 3L], pos, function(t) length(unique(t)))
        if(any(lens > 1L))
            stop("Incompatible constraints.")
        ## Drop duplicated entries.
        A <- A[!ind, , drop = FALSE]
    }
    A
}

### * Incidence generators

### ** .make_incidence_from_upper_tri

.make_incidence_from_upper_tri <-
function(x, family, labels = NULL, diagonal)
{
    ## Compute the indicences of a binary relation from the given family
    ## from its upper triangular part (provided this is possible, of
    ## course).

    family <- match.arg(family, c("E", "L", "S", "T"))

    if(family %in% c("E", "S")) {
        ## Equivalence or symmetric relation.
        y <- diag(diagonal / 2)
        y[upper.tri(y)] <- x
        y <- y + t(y)
    }
    else {
        ## Linear order or tournament.
        y <- diag(diagonal)
        y[upper.tri(y)] <- x
        y[lower.tri(y)] <- 1 - t(y)[lower.tri(y)]
    }

    dimnames(y) <- labels
    y
}

### ** .make_incidence_from_offdiag

.make_incidence_from_offdiag <-
function(x, family, labels = NULL, diagonal)
{
    family <- match.arg(family,
                        .SD_families_using_offdiag_parametrization)

    ## Compute the indicences of a binary relation from its off-diagonal
    ## part (provided this is possible, of course).

    ## <NOTE>
    ## Diagonal entries for the incidences are a mess.
    ## For antisymmetric relations, it really does/should not matter.
    ## We use the standard family definitions to infer reflexivity or
    ## irreflexivity; where neither is implied, we use a majority vote
    ## for reflexivity to determine the diagonal entries.
    ## </NOTE>

    y <- diag(diagonal)
    y[.offdiag(y)] <- x

    dimnames(y) <- labels
    y
}

### ** .make_incidence_from_symdiff_MILP_solution

.make_incidence_from_symdiff_MILP_solution <-
function(x, family, labels = NULL, diagonal)
{
    if(family %in% .SD_families_using_offdiag_parametrization)
        .make_incidence_from_offdiag(round(x),
                                     family, labels, diagonal)
    else
        .make_incidence_from_upper_tri(round(x),
                                       family, labels, diagonal)
}

### ** .make_incidence_from_triples

.make_incidence_from_triples <-
function(x, family, labels = NULL, diagonal)
{
    ## Compute the indicences of a binary relation from a 3-column
    ## matrix the rows of which are triples (i, j, x_{ij}) with x_{ij}
    ## the incidence at position (i, j).

    ## Set up incidences for the non-redundant pairs.
    I <- diag(diagonal)
    I[x[, -3L, drop = FALSE]] <- x[, 3L]
    ## And complete according to family.
    if(family %in% c("E", "L", "S", "T")) {
        ind <- lower.tri(I)
        I[ind] <- if(family %in% c("E", "S"))
            t(I)[ind]
        else
            1 - t(I)[ind]
    }

    dimnames(I) <- labels
    I
}

### ** .make_incidence_from_class_memberships

.make_incidence_from_class_memberships <-
function(M, family, labels)
{
    family <- match.arg(family, c("E", "W"))
    I <- if(family == "E")
        tcrossprod(M)
    else {
        nc <- ncol(M)
        E <- matrix(0, nc, nc)
        E[row(E) <= col(E)] <- 1
        tcrossprod(M %*% E, M)
    }
    dimnames(I) <- labels
    I
}

### ** .make_incidence_from_class_membership_triples

.make_incidence_from_class_membership_triples <-
function(x, family, labels)
{
    M <- matrix(0, max(x[, 1L]), max(x[, 2L]))
    M[x[, -3L, drop = FALSE]] <- x[, 3L]
    ## Could also (more efficiently?) do:
    ##  ind <- which(x[, 3L] > 0)
    ##  M[x[ind, -3L, drop = FALSE]] <- 1
    .make_incidence_from_class_memberships(M, family, labels)
}

### * fit_relation_symdiff_E_k

fit_relation_symdiff_E_k <-
function(C, nc, control = list())
{
    labels <- dimnames(C)

    ## The optimization problem we have to solve is
    ##   \sum_{i,j,k} c_{ij} m_{ik} m_{jk} => max
    ## over all binary stochastic matrices M.
    ## With x = vec(M), this translates into
    ##   x' kronecker(I, C) x => max
    ## under the constraints that x is all binary with 
    ##   kronecker(1', I) x = 1
    ## Rather than simply requiring that all classes are to be used via
    ##   kronecker(I, 1') x >= 1
    ## we actually prefer to arrange them in non-increasing cardinality
    ##   \sum_i m_{i,k} \ge \sum_i m_{i,k+1}, k = 1, ..., nc - 1
    ## and require that the smallest one is non-empty as well:
    ##   \sum_i m_{i,nc} \ge 1

    sparse <- !identical(control$sparse, FALSE)

    no <- nrow(C)
    Q <- kronecker(diag(1, nc, nc), C)
    ## (No need to take halves as linear part is zero.)
    if(sparse)
        Q <- as.simple_triplet_matrix(Q)
    mat <- rbind(kronecker(rbind(rep.int(1, nc)), diag(1, no, no)),
                 ## The first matrix has 1 in the main and -1 in the
                 ## first upper diagonal.
                 kronecker(matrix(c(rep.int(c(1,
                                              rep.int(0, nc - 1L),
                                              -1),
                                            nc - 1L),
                                    1),
                                  nrow = nc, ncol = nc),
                           rbind(rep.int(1, no))))
    ## (Of course, we can generate mat more efficiently ...)
    if(sparse)
        mat <- as.simple_triplet_matrix(mat)
    dir <- rep.int(c("==", ">="), c(no, nc))
    rhs <- rep.int(c(1, 0, 1), c(no, nc - 1L, 1L))

    ## Handle possibly additional explicit constrains.
    ## (Which specify that pairs of objects are in relation or not.)
    if(!is.null(A <- control$constraints)) {
        pos <- function(i, k) {
            ## Position of variable m_{ik}.
            i + (k - 1L) * no
        }
        acmaker <-
            .make_additional_constraint_maker_using_memberships_E(pos,
                                                                  no * nc,
                                                                  nc, sparse)
        A <- .canonicalize_additional_constraints(A, "E")
        add <- acmaker(A)
        mat <- rbind(mat, add$mat)
        dir <- c(dir, add$dir)
        rhs <- c(rhs, add$rhs)
    }

    nos <- .n_from_control_list(control)
    one <- nos == 1L
    ## Membership matrices of equivalence relations are only unique up
    ## to column permutations.  At worst, all classes have the same
    ## number of elements, so arranging in non-increasing cardinality
    ## does not reduce the number of solutions found.
    nom <- as.integer(min(.Machine$integer.max, nos * factorial(nc)))
    solver <- control$solver
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")
    ## Empirical evidence suggests that explicit linearization works
    ## faster even for cplex, so linearize by default.
    linearize <- !identical(control$linearize, FALSE)
    control <- c(list(n = nom, verbose = verbose, linearize = linearize),
                 control$control)
    if(!one) {
        ## Split row-wise (as the first 1 in a row implies all other
        ## entries are 0).
        order <- c(matrix(seq_len(no * nc), nrow = nc, ncol = no,
                          byrow = TRUE))
        control$order <- order
    }
    out <- solve_MIQP(MIQP(list(Q, double(nrow(Q))),
                           list(mat = mat, dir = dir, rhs = rhs),
                           types = "B", maximum = TRUE),
                      solver, control)
    finisher <- function(e) {
        .stop_if_lp_status_is_nonzero(e$status, "E")
        M <- matrix(e$solution, ncol = nc)
        .make_incidence_from_class_memberships(M, "E", labels)
    }
    if(!one) {
        ## Note that equivalence classes are only unique up to
        ## permutations of the class labels.
        out <- unique(lapply(out, finisher))
        if(length(out) > nos)
            out <- out[seq_len(nos)]
        out
    }
    else finisher(out)

}

### * fit_relation_symdiff_W_k

fit_relation_symdiff_W_k <-
function(C, nc, control = list())
{
    labels <- dimnames(C)

    ## The optimization problem we have to solve is
    ##   \sum_{i,j,k,l} c_{ij} I(k <= l) m_{ik} m_{jl} => max
    ## over all binary stochastic matrices M.
    ## With x = vec(M), this translates into
    ##   x' kronecker(J, C) x => max
    ## (where J_{kl} is one iff k <= l) under the constraints that x is
    ## all binary with 
    ##   kronecker(1', I) x = 1
    ## and if all classes are to be used,
    ##   kronecker(I, 1') x >= 1

    sparse <- !identical(control$sparse, FALSE)

    J <- matrix(0, nc, nc)
    J[row(J) <= col(J)] <- 1
    no <- nrow(C)
    Q <- kronecker(J, C)
    ## (No need to take halves as linear part is zero.)
    if(sparse)
        Q <- as.simple_triplet_matrix(Q)
    mat <- rbind(kronecker(rbind(rep.int(1, nc)), diag(1, no, no)),
                 kronecker(diag(1, nc, nc), rbind(rep.int(1, no))))
    ## (Of course, we can generate mat more efficiently ...)
    if(sparse)
        mat <- as.simple_triplet_matrix(mat)
    dir <- rep.int(c("==", ">="), c(no, nc))
    rhs <- rep.int(1, nc + no)

    ## Handle constraints on the class sizes.
    if(!is.null(l <- control$l)) {
        ## Verify that l is feasible.
        ## It would be nice to be nice, in particular by rescaling the
        ## sizes to sum to the number of objects (so that one can
        ## e.g. specify proportions).  But using something like
        ##   l <- round(l / sum(l) * nc)
        ## does not work (e.g., c(2.4, 2.4, 5.2) sums to 10 but after
        ## rounding only to 9).
        if((length(l) != nc) || (sum(l) != no))
            stop("Invalid specification of class sizes.")
        ## And now add constraints \sum_i m_{ik} = l_k, k = 1, ..., nc.
        mat <- rbind(mat,
                     kronecker(diag(1, nc), rbind(rep.int(1, no))))
        dir <- c(dir, rep.int("==", nc))
        rhs <- c(rhs, l)
    }
    ## Handle additional explicit constrains.
    ## (Which specify that pairs of objects are in relation or not.)
    if(!is.null(A <- control$constraints)) {
        pos <- function(i, k) {
            ## Position of variable m_{ik}.
            i + (k - 1L) * no
        }
        acmaker <-
            .make_additional_constraint_maker_using_memberships_W(pos,
                                                                  no * nc,
                                                                  nc, sparse)
        A <- .canonicalize_additional_constraints(A, "W")
        add <- acmaker(A)
        mat <- rbind(mat, add$mat)
        dir <- c(dir, add$dir)
        rhs <- c(rhs, add$rhs)
    }

    nos <- .n_from_control_list(control)
    one <- nos == 1L
    solver <- control$solver
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")
    ## Empirical evidence suggests that explicit linearization works
    ## faster even for cplex, so linearize by default.
    linearize <- !identical(control$linearize, FALSE)
    control <- c(list(n = nos, verbose = verbose, linearize = linearize),
                 control$control)
    if(!one) {
        ## Split row-wise (as the first 1 in a row implies all other
        ## entries are 0).
        order <- c(matrix(seq_len(no * nc), nrow = nc, ncol = no,
                          byrow = TRUE))
        control$order <- order
    }

    out <- solve_MIQP(MIQP(list(Q, double(nrow(Q))),
                           list(mat = mat, dir = dir, rhs = rhs),
                           types = "B", maximum = TRUE),
                      solver, control)
    finisher <- function(e) {
        .stop_if_lp_status_is_nonzero(e$status, "W")
        M <- matrix(e$solution, ncol = nc)
        .make_incidence_from_class_memberships(M, "W", labels)
    }
    if(!one) lapply(out, finisher) else finisher(out)
}

### * fit_relation_CKS_via_QP

fit_relation_CKS_via_QP <-
function(P, Q, family, control)
{
    ## Note: relevant families are those which are neither complete nor
    ## symmetric nor antisymmetric (currently, only transitive relations
    ## and preorders).

    ## Set up coefficients of binary quadratic program.
    A <- Q - P - t(P)
    A[row(A) >= col(A)] <- 0
    B <- P - Q
    diag(B) <- - diag(Q)

    ## It would be "natural" to rewrite the objective function in terms
    ## of vec([x_{ij}]).  But we have transitivity constraint makers
    ## only for offdiag (and upper.tri) parametrizations, so we use the
    ## former and deal directly with the diagonal terms.

    sparse <- !identical(control$sparse, FALSE)

    n <- nrow(B)
    NP <- .n_of_pairs(n, family)
    pos <- .make_pos(n, family)
    mat <- .make_transitivity_constraint_mat(n, family, sparse)
    dir <- .make_transitivity_constraint_dir(n, family)
    rhs <- .make_transitivity_constraint_rhs(n, family)
    ## As we also allow for the families using offdiag parametrizatios
    ## for which the consensus problem admits an LP formulation to use
    ## the QP formulation: 
    if(! family %in% .SD_families_using_offdiag_parametrization) {
        ## Unfortunately, QP formulations for families using upper_tri
        ## paramatrization are not supported.  Doable via maintaining
        ## the offdiag parametrization and adding constraints for
        ## symmetry or antisymmetry and generating full transitivity
        ## constraints where needed ...
        stop("Not implemented.")
    } else if(family %in% .SD_families_using_tot_or_asy_constraints) {
        mat <- rbind(mat, .make_tot_or_asy_constraint_mat(n, sparse))
        dir <- c(dir, .make_tot_or_asy_constraint_dir(n, family))
        rhs <- c(rhs, .make_tot_or_asy_constraint_rhs(n))
    }
                     
    ## When using the QP formulation, need to figure out how to rewrite
    ##  \sum_{i,j: i < j} \alpha_{ij} x_{ij} x_{ji}
    ##    + \sum_{ij} \beta_{ij} x_{ij}.
    ## in terms of vec(offdiag([x_{ij}])).  Note that using a dense QP
    ## formulation is not very efficient as we really only have O(n^2)
    ## quadratic terms (but a dense Q has O(n^4) elements).
    ind <- which(row(A) < col(A), arr.ind = TRUE)
    i <- ind[, 1L]
    j <- ind[, 2L]
    if(sparse) {
        Q <- simple_triplet_matrix(pos(i, j), pos(j, i), A[ind],
                                   nrow = NP, ncol = NP)
    } else {
        Q <- matrix(0, NP, NP)
        Q[cbind(pos(i, j), pos(j, i))] <- A[ind]
    }

    nos <- .n_from_control_list(control)
    one <- nos == 1L
    solver <- control$solver
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")
    ## Empirical evidence suggests that explicit linearization works
    ## faster even for cplex, so linearize by default.
    linearize <- !identical(control$linearize, FALSE)
    control <- c(list(n = nos, verbose = verbose, linearize = linearize),
                 control$control)
    
    out <- solve_MIQP(MIQP(list(2 * Q, B[.offdiag(B)]),
                           list(mat = mat, dir = dir, rhs = rhs),
                           types = "B", maximum = TRUE),
                      solver, control)
    ## Simplify matters by *not* expanding all diagonal solutions.
    diagonal <- if(family == "preorder")
        rep.int(1, n)
    else
        diag(B) >= 0
    finisher <- function(e) {
        if(e$status != 0)
            stop("MIQP could not be solved.")
        .make_incidence_from_offdiag(round(e$solution),
                                     family, dimnames(P), diagonal)
    }
    if(!one) lapply(out, finisher) else finisher(out)
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
