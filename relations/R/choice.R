relation_choice <-
function(x, method = "symdiff", weights = 1, control = list(), ...)
{
    dots <- list(...)
    control[names(dots)] <- dots

    if(inherits(x, "gset")) {
        relations <- relation_ensemble(list = gset_support(x))
        control$weights <- gset_memberships(x)
    } else {
        relations <- as.relation_ensemble(x)
    }

    if(!length(relations))
        stop("Cannot compute choice from empty ensemble.")

    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")
    
    weights <- rep(weights, length.out = length(relations))
    if(any(weights < 0))
        stop("Argument 'weights' has negative elements.")
    if(!any(weights > 0))
        stop("Argument 'weights' has no positive elements.")

    known_methods <-
        list("symdiff" = ".relation_choice_symdiff",
             "euclidean" = ".relation_choice_euclidean",
             "CKS" = ".relation_choice_CKS",
             "Schulze" = ".relation_choice_Schulze"
             )
    if(is.character(method)) {
        ## Hopefully of length one, add some tests eventually ...
        if(is.na(ind <- pmatch(method, names(known_methods))))
            stop(gettextf("Method '%s' is not a valid choice method."),
                 domain = NA)
        method <- get(known_methods[[ind]][1L])
    }
    else if(!is.function(method))
        stop("Invalid 'method' argument.")

    method(relations, weights, control)
}

## <FIXME>
## Add information on which problem is solved, and how this is done.
## </FIXME>

.relation_choice_symdiff <-
function(relations, weights, control, euclidean = FALSE)
{
    if(!euclidean && !.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    ## Argument handling.
    k <- control$k
    if(is.null(k)) k <- 1L              # Single winner by default.

    ## Coefficients of the choice QP.
    B <- .make_fit_relation_symdiff_M(relations, weights)
    if(identical(control$reverse, TRUE))
        B <- t(B)
    ## Explicitly, we have
    ##   C <- pmax(B, 0)
    ##   M <- C + t(C) - B
    ##   diag(M) <- 0
    ## which can be simplified to:    
    M <- pmax(t(B), 0) - pmin(B, 0)
    diag(M) <- 0

    ## Underlying set of objects to choose from.
    ## (Domain of the choice problem.)
    D <- .get_elements_in_homorelation(relations)

    .find_SD_or_CKS_choice(M, k, D, control)
}

.relation_choice_euclidean <-
function(relations, weights, control)
    .relation_choice_symdiff(relations, weights, control, TRUE)

.relation_choice_CKS <-
function(relations, weights, control)    
{
    if(!.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    ## Argument handling.
    k <- control$k
    if(is.null(k)) k <- 1L              # Single winner by default.

    ## Coefficients of the choice QP.
    incidences <- lapply(relations, relation_incidence)
    P <- .make_fit_relation_CKS_P(incidences, weights)
    if(identical(control$reverse, TRUE))
        P <- t(P)
    Q <- .make_fit_relation_CKS_Q(incidences, weights)
    ## Explicitly, we have
    ##   B <- P - Q
    ##   diag(B) <- - diag(Q)
    ##   C <- (pmax(Q, P, t(P), 0) - Q) / 2
    ##   diag(C) <- pmax(- diag(Q), 0)
    ##   M <- 2 * C - B
    ##   diag(M) <- 0
    ## which can be simplified to:
    M <- pmax(0, P, t(P), Q) - P
    diag(M) <- 0

    ## Underlying set of objects to choose from.
    ## (Domain of the choice problem.)
    D <- .get_elements_in_homorelation(relations)

    .find_SD_or_CKS_choice(M, k, D, control)
}

.find_SD_or_CKS_choice <-
function(M, k, D, control)
{
    ## It is somewhat unclear which formulation in general works
    ## better.
    ## For CPLEX, it seems that any linearization of the BQP is better
    ## than directly solving the QP formulation.
    ## For the other solvers, the "tailor-made" direct linearization
    ## seems to perform better than the general purpose linearization of
    ## the QP formulation.  Not sure why, as in general the latter
    ## results in smaller problems---its number of continuous variables
    ## is n(n-1)/2 minus half the number of zero off-diagonal terms of
    ## sym(M), whereas for the former the number is n(n-1) minus the
    ## number of zero off-diagonal terms.
    ## But for the direct linearization the binary variables occur in
    ## the constraints only---maybe this accounts for the difference.
    ## Hence, by default we use the direct linearization.
    ## We might change the QP formulation to also linearize by default.
    if(identical(control$QP, TRUE))
        .find_SD_or_CKS_choice_via_QP(M, k, D, control)
    else
        .find_SD_or_CKS_choice_via_LP(M, k, D, control)
}

## Solve SD/CKS choice problem via LP.

.find_SD_or_CKS_choice_via_LP <-
function(M, k, D, control)
{
    ## Solve SD/CKS choice problem via LP using a "tailor-made"
    ## linearization.
    ## With u the committee/winners indicator, the BQP is
    ##   \sum_{i,j} (1 - u_i) m_{ij} u_j \to \min
    ## subject to \sum_i u_i = k, where m_{ij} \ge 0 and m_{ii} = 0.
    ## Letting I be the set of (i,j) with m_{ij} > 0 and
    ##   z_{ij} = (1 - u_i) u_j
    ## for (i,j) \in I, we can write the BQP as
    ##   \sum_I m_{ij} z_{ij} \to \min
    ## subject to
    ##   z_{ij} + u_i - u_j \ge 0
    ## (as z_{ij} \ge (1 - u_i) + u_j - 1 = - u_i + u_j) and
    ##   \sum_i u_i = k.
    ## We use c(u, z) for the decision variables.

    ## Argument handling.
    nos <- .n_from_control_list(control)
    MIP <- identical(control$MIP, TRUE)
    sparse <- !identical(control$sparse, FALSE)
    solver <- control$solver
    control <- control$control

    ind <- which(M > 0, arr.ind = TRUE)
    n_u <- nrow(M)
    n_z <- nrow(ind)
    pos <- seq_len(n_z)
    if(sparse) {
        mat <- simple_triplet_matrix(c(rep.int(pos, 3L),
                                       rep.int(n_z + 1L, n_u)),
                                     c(c(ind), n_u + pos, seq_len(n_u)),
                                     rep.int(c(1, -1, 1),
                                             c(n_z, n_z, n_z + n_u)),
                                     n_z + 1L,
                                     n_u + n_z)
    } else {
        mat <- matrix(0, n_z, n_u)
        mat[cbind(pos, ind[, 1L])] <- 1
        mat[cbind(pos, ind[, 2L])] <- -1
        mat <- rbind(cbind(mat, diag(1, n_z)),
                     rep.int(c(1, 0), c(n_u, n_z)))
    }

    milp <- MILP(c(rep.int(0, n_u), M[ind]),
                 list(mat,
                      c(rep.int(">=", n_z), "=="),
                      c(double(n_z), k)),
                 types = rep.int(c("B", "C"), c(n_u, n_z)),
                 maximum = FALSE)
    out <- solve_MILP(milp, solver, c(list(n = nos), control))
    
    pos <- seq_len(n_u)  
    finisher <- function(e) as.set(D[e$solution[pos] == 1])
    out <- if(nos > 1L) lapply(out, finisher) else finisher(out)
    if(MIP) attr(out, "MIP") <- milp
    out
}

## Solve SD/CKS choice problem via QP.

.find_SD_or_CKS_choice_via_QP <-
function(M, k, D, control)
{
    ## Argument handling.
    nos <- .n_from_control_list(control)
    MIP <- identical(control$MIP, TRUE)
    sparse <- !identical(control$sparse, FALSE)
    solver <- control$solver
    control <- control$control

    n <- nrow(M)
    mat <- if(sparse)
        simple_triplet_matrix(rep.int(1, n), seq_len(n), rep.int(1, n),
                              1L, n)
    else
        matrix(1, 1L, n)

    miqp <- MIQP(list(2 * M, - colSums(M)),
                 list(mat, "==", k),
                 types = "B",
                 maximum = TRUE)
    out <- solve_MIQP(miqp, solver, c(list(n = nos), control))
    
    finisher <- function(e) as.set(D[e$solution == 1])
    out <- if(nos > 1L) lapply(out, finisher) else finisher(out)
    if(MIP) attr(out, "MIP") <- miqp
    out
}

.relation_choice_Schulze <-
function(relations, weights, control)
{
    if(!.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    ## Need the numbers of voters who strictly prefer candidate i to
    ## candidate j (i > j, i.e. not(i <= j)): we compute B - d

    d <- .weighted_sum_of_arrays(lapply(relations, .incidence), weights)

    ## From which we compute
    ##   p[i,j] = d[i,j] - d[j,i]
    p <- t(d) - d

    N <- NROW(p)

    ## The reference has
    ## for(i in seq_len(N)) {
    ##     for(j in seq_len(N)) {
    ##         if(i != j) {
    ##             for(k in seq_len(N)) {
    ##                 if((i != k) && (j != k)) {
    ##                     s <- min(p[j, i], p[i, k])
    ##                     if(p[j, k] < s)
    ##                         p[j, k] <- s
    ##                 }
    ##             }
    ##         }
    ##     }
    ## }
    ## but it seems we can simply ignore the inequalities ...

    for(i in seq_len(N))
        p <- pmax(p, outer(p[, i], p[i, ], pmin))

    ## Determine potential winners.

    D <- .get_elements_in_homorelation(relations)
    
    as.set(D[apply(p >= t(p), 1L, all)])

    ## Or equivalently:
    ##   ! apply(t(p) > p, 1L, any)
}
