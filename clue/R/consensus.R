### * cl_consensus

cl_consensus <-
function(x, method = NULL, weights = 1, control = list())
{
    ## <NOTE>
    ## Interfaces are a matter of taste.
    ## E.g., one might want to have a 'type' argument indication whether
    ## hard or soft partitions are sought.  One could then do
    ##   cl_consensus(x, method = "euclidean", type = "hard")
    ## to look for an optimal median (or least squares) hard partition
    ## (for euclidean dissimilarity).
    ## For us, "method" really indicates a certain algorithm, with its
    ## bells and whistles accessed via the 'control' argument.
    ## </NOTE>

    clusterings <- as.cl_ensemble(x)

    if(!length(clusterings))
        stop("Cannot compute consensus of empty ensemble.")

    weights <- rep(weights, length.out = length(clusterings))
    if(any(weights < 0))
        stop("Argument 'weights' has negative elements.")
    if(!any(weights > 0))
        stop("Argument 'weights' has no positive elements.")

    if(!is.function(method)) {
        if(!inherits(method, "cl_consensus_method")) {
            ## Get the method definition from the registry.
            type <- .cl_ensemble_type(clusterings)
            if(is.null(method))
                method <- .cl_consensus_method_default(type)                
            method <- get_cl_consensus_method(method, type)            
        }
        method <- method$definition
    }

    method(clusterings, weights, control)
}

### * .cl_consensus_partition_DWH

.cl_consensus_partition_DWH <-
function(clusterings, weights, control)
{
    ## <TODO>
    ## Could make things more efficient by subscripting on positive
    ## weights.
    ## (Note that this means control$order has to be subscripted as
    ## well.)
    ## </TODO>
    
    max_n_of_classes <- max(sapply(clusterings, n_of_classes))

    ## Control parameters.
    k <- control$k
    if(is.null(k))
        k <- max_n_of_classes
    order <- control$order
    if(is.null(order))
        order <- sample(seq_along(clusterings))

    clusterings <- clusterings[order]
    weights <- weights[order]

    k_max <- max(k, max_n_of_classes)
    s <- weights / cumsum(weights)
    s[is.na(s)] <- 0                    # Division by zero ...
    
    M <- cl_membership(clusterings[[1L]], k_max)
    for(b in seq_along(clusterings)[-1L]) {
        mem <- cl_membership(clusterings[[b]], k_max)
        ## Match classes from conforming memberships.
        ind <- solve_LSAP(crossprod(M, mem), maximum = TRUE)
        M <- (1 - s[b]) * M + s[b] * mem[, ind]
        if(k < k_max)
            M <- .project_to_leading_columns(M, k)
    }

    M <- .cl_membership_from_memberships(M[, seq_len(k), drop = FALSE], k)

    as.cl_partition(M)
}

### * .cl_consensus_partition_AOS

.cl_consensus_partition_AOS <-
function(clusterings, weights, control,
         type = c("SE", "HE", "SM", "HM"))
{
    ## The start of a general purpose optimizer for determining
    ## consensus partitions by minimizing
    ##   \sum_b w_b d(M, M_b) ^ e
    ##     = \sum_b \min_{P_b} w_b f(M, M_b P_b) ^ e
    ## for the special case where the criterion function is based on
    ## M and M_b P_b (i.e., column permutations of M_b), as opposed to
    ## the general case where d(M, M_b) = \min_{P_b} f(M, P_b, M_b)
    ## handled by .cl_consensus_partition_AOG().
    ##
    ## The AO ("alternative optimization") proceeds by alternatively
    ## matching the M_b to M by minimizing f(M, M_b P_b) over P_b, and
    ## fitting M by minimizing \sum_b w_b f(M, M_b P_b) ^ e for fixed
    ## matchings.
    ##
    ## Such a procedure requires three ingredients: a function for
    ## matching M_b to M (in fact simply replacing M_b by the matched
    ## M_b P_b); a function for fitting M to the \{M_b P_b\}, and a
    ## function for computing the value of the criterion function
    ## corresponding to this fit (so that one can stop if the relative
    ## improvement is small enough).
    ##
    ## For the time being, we only use this to determine soft and hard
    ## Euclidean least squares consensus partitions (soft and hard
    ## Euclidean means), so the interface does not yet reflect the
    ## generality of the approach (which would either pass the three
    ## functions, or even set up family objects encapsulating the three
    ## functions).
    ##
    ## This special case is provided for efficiency and convenience.
    ## Using the special form of the criterion function, we can simply
    ## always work memberships with the same maximal number of columns,
    ## and with the permuted \{ M_b P_b \}.

    ## For the time being ...
    type <- match.arg(type)

    w <- weights / sum(weights)    
    n <- n_of_objects(clusterings)
    k_max <- max(sapply(clusterings, n_of_classes))

    ## Control parameters.
    k <- control$k
    if(is.null(k))
        k <- k_max
    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 100
    nruns <- control$nruns
    reltol <- control$reltol
    if(is.null(reltol))
        reltol <- sqrt(.Machine$double.eps)
    start <- control$start
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")

    ## Handle start values and number of runs.
    if(!is.null(start)) {
        if(!is.list(start)) {
            ## Be nice to users.
            start <- list(start)
        }
        nruns <- length(start)
    } else {
        if(is.null(nruns)) {
            ## Use nruns only if start is not given.
            nruns <- 1L
        }
        start <- replicate(nruns,
                           .random_stochastic_matrix(n, k),
                           simplify = FALSE)
    }

    ## The maximal (possible) number of classes in M and the \{ M_b \}.
    k_all <- max(k, k_max)

    value <-
        switch(type,
               SE = , HE = function(M, memberships, w) {
                   sum(w * sapply(memberships,
                                  function(u) sum((u - M) ^ 2)))
               },
               SM = , HM = function(M, memberships, w) {
                   sum(w * sapply(memberships,
                                  function(u) sum(abs(u - M))))
               })
    ## Return the M[, ind] column permutation of M optimally matching N.
    match_memberships <-
        switch(type,
               SE = , HE = function(M, N) {
                   M[, solve_LSAP(crossprod(N, M), maximum = TRUE),
                     drop = FALSE]
               },
               SM = , HM = function(M, N) {
                   M[, solve_LSAP(.cxdist(N, M, "manhattan")),
                     drop = FALSE]
               })
    ## Function for fitting M to (fixed) memberships \{ M_b P_b \}.
    ## As we use a common number of columns for all membership matrices
    ## involved, we need to pass the desired 'k' ...
    fit_M <-
        switch(type,
               SE = function(memberships, w, k) {
                   ## Update M as \sum w_b M_b P_b.
                   M <- .weighted_sum_of_matrices(memberships, w, nrow(M))
                   ## If k < k_all, "project" as indicated in Gordon &
                   ## Vichi (2001), p. 238.
                   if(k < ncol(M))
                       M <- .project_to_leading_columns(M, k)
                   M
               },
               HE = , HM = function(memberships, w, k) {
                   ## Compute M as \sum w_b M_b P_b.
                   M <- .weighted_sum_of_matrices(memberships, w, nrow(M))
                   ## And compute a closest hard partition H(M) from
                   ## that, using the first k columns of M.
                   ids <- max.col(M[ , seq_len(k), drop = FALSE])
                   .cl_membership_from_class_ids(ids, ncol(M))
               },
               SM = .l1_fit_M)

    memberships <- lapply(clusterings, cl_membership, k_all)

    V_opt <- Inf
    M_opt <- NULL
    for(run in seq_along(start)) {
        if(verbose && (nruns > 1L))
            message(gettextf("AOS run: %d", run))
        M <- start[[run]]
        if(k < k_all)
            M <- cbind(M, matrix(0, nrow(M), k_all - k))
        memberships <- lapply(memberships, match_memberships, M)
        old_value <- value(M, memberships, w)
        if(verbose)
            message(gettextf("Iteration: 0 *** value: %g", old_value))
        iter <- 1L
        while(iter <= maxiter) {
            ## Fit M to the M_b P_b.
            M <- fit_M(memberships, w, k)
            ## Match the \{ M_b P_b \} to M.
            memberships <- lapply(memberships, match_memberships, M)
            ## Update value.
            new_value <- value(M, memberships, w)
            if(verbose)
                message(gettextf("Iteration: %d *** value: %g",
                                 iter, new_value))
            if(abs(old_value - new_value)
               < reltol * (abs(old_value) + reltol))
                break
            old_value <- new_value
            iter <- iter + 1L
        }
        if(new_value < V_opt) {
            converged <- (iter <= maxiter)
            V_opt <- new_value
            M_opt <- M
        }
        if(verbose)
            message(gettextf("Minimum: %g", V_opt))
    }

    M <- .stochastify(M_opt)
    rownames(M) <- rownames(memberships[[1L]])
    meta <- list(objval = value(M, memberships, w),
                 converged = converged)
    M <- .cl_membership_from_memberships(M[, seq_len(k), drop = FALSE],
                                         k, meta)
    
    as.cl_partition(M)
}

.random_stochastic_matrix <-
function(n, k)
{
    M <- matrix(runif(n * k), n, k)
    M / rowSums(M)
}

.l1_fit_M <-
function(memberships, w, k)
{
    ## Determine stochastic matrix M with at most k leading nonzero
    ## columns such that
    ##
    ##    \sum_b w_b \sum_{i,j} | m_{ij}(b) - m_{ij} | => min
    ##
    ## where the sum over j goes from 1 to k.
    ##
    ## Clearly, this can be done separately for each row, where we need
    ## to minimize
    ##
    ##    \sum_b w_b \sum_j | y_j(b) - x_j | => min
    ##
    ## over all probability vectors x.  Such problems can e.g. be solved
    ## via the following linear program:
    ##
    ##    \sum_b \sum_j w_b e'(u(b) + v(b)) => min
    ##
    ## subject to
    ##
    ##    u(1), v(1), ..., u(B), v(B), x >= 0
    ##                   x + u(b) - v(b)  = y(b),    b = 1, ..., B
    ##                               e'x  = 1
    ##
    ## (where e = [1, ..., 1]).
    ##
    ## So we have one long vector z of "variables":
    ##
    ##    z = [u(1)', v(1)', ..., u(B)', v(B)', x']'
    ##
    ## of length (2B + 1) k, with x the object of interest.
    
    ## Rather than providing a separate function for weighted L1 fitting
    ## of probability vectors we prefer doing "everything" at once, in
    ## order to avoid recomputing the coefficients and constraints of
    ## the associated linear program.

    B <- length(memberships)
    L <- (2 * B + 1) * k

    ## Set up associated linear program.

    ## Coefficients in the objective function.
    objective_in <- c(rep(w, each = 2 * k), rep.int(0, k))
    
    ## Constraints.
    constr_mat <- rbind(diag(1, L),
                        cbind(kronecker(diag(1, B),
                                        cbind(diag(1, k),
                                              diag(-1, k))),
                              kronecker(rep.int(1, B),
                                        diag(1, k))),
                        c(rep.int(0, 2 * B * k), rep.int(1, k)))
    constr_dir <- c(rep.int(">=", L), rep.int("==", B * k + 1L))

    ind <- seq.int(from = 2 * B * k + 1L, length.out = k)
    nr <- NROW(memberships[[1L]])
    nc <- NCOL(memberships[[1L]])
    M <- matrix(0, nrow = nr, ncol = k)

    ## Put the memberships into one big array so that we can get their
    ## rows more conveniently (and efficiently):

    memberships <- array(unlist(memberships), c(nr, nc, B))

    for(i in seq_len(nr)) {
        out <- lpSolve::lp("min",
                           objective_in,
                           constr_mat,
                           constr_dir,
                           c(rep.int(0, L), memberships[i, seq_len(k), ], 1))
        M[i, ] <- out$solution[ind]
    }

    ## Add zero columns if necessary.
    if(k < nc)
        M <- cbind(M, matrix(0, nr, nc - k))

    M
}

### ** .cl_consensus_partition_soft_euclidean

.cl_consensus_partition_soft_euclidean <-
function(clusterings, weights, control)
    .cl_consensus_partition_AOS(clusterings, weights, control, "SE")

### ** .cl_consensus_partition_hard_euclidean

.cl_consensus_partition_hard_euclidean <-
function(clusterings, weights, control)
    .cl_consensus_partition_AOS(clusterings, weights, control, "HE")

### ** .cl_consensus_partition_soft_manhattan

.cl_consensus_partition_soft_manhattan <-
function(clusterings, weights, control)
    .cl_consensus_partition_AOS(clusterings, weights, control, "SM")

### ** .cl_consensus_partition_hard_manhattan

.cl_consensus_partition_hard_manhattan <-
function(clusterings, weights, control)
    .cl_consensus_partition_AOS(clusterings, weights, control, "HM")

### * .cl_consensus_partition_AOG

.cl_consensus_partition_AOG <-
function(clusterings, weights, control, type = c("GV1"))
{
    ## The start of a general purpose optimizer for determining
    ## consensus partitions by minimizing
    ##   \sum_b w_b d(M, M_b) ^ p
    ##     = \sum_b \min_{P_b} w_b f(M, M_b, P_b) ^ e
    ## for general dissimilarity matrices which involve class matching
    ## via permutation matrices P_b.
    ##
    ## The AO ("Alternative Optimization") proceeds by alternating
    ## between determining the optimal permutations P_b by minimizing
    ##   f(M, M_b, P_b)
    ## for fixed M, and fitting M by minimizing
    ##   \sum_b w_b f(M, M_b, P_b) ^ e
    ## for fixed \{ P_b \}.
    ##
    ## We encapsulate this into functions fit_P() and fit_M() (and a
    ## value() function for the criterion function to be minimized with
    ## respect to both M and \{ P_b \}, even though the current
    ## interface does not yet reflect the generality of the approach.
    ##
    ## Note that rather than passing on information about the numbers of
    ## classes (e.g., needed for GV1) and representing all involved
    ## membership matrices with the same maximal number of columns, we
    ## use "minimal" representations with no dummy classes (strictly
    ## speaking, with the possible exception of M, for which the given k
    ## is used).

    ## For the time being ...
    type <- match.arg(type)

    w <- weights / sum(weights)    
    n <- n_of_objects(clusterings)
    k_max <- max(sapply(clusterings, n_of_classes))
    
    ## Control parameters.
    k <- control$k
    if(is.null(k))
        k <- k_max
    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 100L
    nruns <- control$nruns
    reltol <- control$reltol
    if(is.null(reltol))
        reltol <- sqrt(.Machine$double.eps)
    start <- control$start
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")

    ## Handle start values and number of runs.
    if(!is.null(start)) {
        if(!is.list(start)) {
            ## Be nice to users.
            start <- list(start)
        }
        nruns <- length(start)
    } else {
        if(is.null(nruns)) {
            ## Use nruns only if start is not given.
            nruns <- 1L
        }
        start <- replicate(nruns,
                           .random_stochastic_matrix(n, k),
                           simplify = FALSE)
    }

    ## <NOTE>
    ## For the given memberships, we can simply use ncol() in the
    ## computations (rather than n_of_classes(), because we used
    ## cl_membership() to create them.  For M, the number of classes
    ## could be smaller than the given k "target".
    ## </NOTE>
    
    value <- function(M, permutations, memberships, w) {

        k <- .n_of_nonzero_columns(M)
        d <- function(u, p) {
            ## Compute the squared GV1 dissimilarity between M and u
            ## based on the M->u class matching p.
            nc_u <- ncol(u)
            if(nc_u == k) {
                ## Simple case: all classes are matched.
                sum((u[, p] - M) ^ 2)
            }
            else {
                ## Only include the matched non-dummy classes of M ..
                ind <- seq_len(k)
                ## ... which are matched to non-dummy classes of u.
                ind <- ind[p[ind] <= nc_u]
                sum((u[, p[ind]] - M[, ind]) ^ 2)
            }
        }
        
        sum(w * mapply(d, memberships, permutations))
    }

    fit_P <- function(u, M) {
        
        ## Return a permutation representing a GV1 optimal matching of
        ## the columns of M to the columns of u (note the order of the
        ## arguments), using a minimal number of dummy classes (i.e., p
        ## has max(.n_of_nonzero_columns(M), n_of_classes(u)) entries).

        ## See also .cl_dissimilarity_partition_GV1().
        
        C <- outer(colSums(M ^ 2), colSums(u ^ 2), "+") -
            2 * crossprod(M, u)
        nc_M <- .n_of_nonzero_columns(M)
        nc_u <- ncol(u)
        ## (See above for ncol() vs n_of_classes().)
        if(nc_M < nc_u)
            C <- rbind(C, matrix(0, nrow = nc_u - nc_M, ncol = nc_u))
        else if(nc_M > nc_u)
            C <- cbind(C, matrix(0, nrow = nc_M, ncol = nc_M - nc_u))
        
        solve_LSAP(C)
    }

    fit_M <- function(permutations, memberships, w) {
        
        ## Here comes the trickiest part ...
        ##
        ## In general, M = [m_{iq}] is determined as follows.
        ## Write value(M, permutations, memberships, w) as
        ##   \sum_b \sum_i \sum_{p=1}^{k_b} \sum_{q=1}^k
        ##      w_b (u_{ip}(b) - m_{iq})^2 x_{pq}(b)
        ## where U(b) and X(b) are the b-th membership matrix and the
        ## permutation matrix representing the M->U(b) non-dummy class
        ## matching (as always, note the order of the arguments).
        ##
        ## Let
        ##   \beta_{iq} = \sum_b \sum_{p=1}^{k_b} w_b u_{ip}(b) x_{pq}(b)
        ##   \alpha_q   = \sum_b \sum_{p=1}^{k_b} w_b x_{pq}(b)
        ## and
        ##   \bar{m}_{iq} =
        ##     \cases{\beta_{iq}/\alpha_q, & $\alpha_q > 0$ \cr
        ##            0                    & otherwise}.
        ## Then, as the cross-product terms cancel out, the value
        ## function rewrites as
        ##   \sum_b \sum_i \sum_{p=1}^{k_b} \sum_{q=1}^k
        ##      w_b (u_{ip}(b) - \bar{m}_{iq})^2 x_{pq}(b)
        ##   + \sum_i \sum_q \alpha_q (\bar{m}_{iq} - m_{iq}) ^ 2,
        ## where the first term is a constant, and the minimum is found
        ## by solving
        ##   \sum_q \alpha_q (\bar{m}_{iq} - m_{iq}) ^ 2 => min!
        ## s.t.
        ##   m_{i1}, ..., m_{ik} >= 0, \sum_{iq} m_{iq} = 1.
        ##
        ## We can distinguish three cases.
        ## A. If S_i = \sum_q \bar{m}_{iq} = 1, things are trivial.
        ## B. If S_i = \sum_q \bar{m}_{iq} < 1.
        ##    B1. If some \alpha_q are zero, then we can choose
        ##        m_{iq} = \bar{m}_{iq} for those q with \alpha_q = 0;
        ##        m_{iq} = 1 / number of zero \alpha's, otherwise.
        ##    B2. If all \alpha_q are positive, we can simply
        ##        equidistribute 1 - S_i over all classes as written
        ##        in G&V.
        ## C. If S_i > 1, things are not so clear (as equidistributing
        ##    will typically result in violations of the non-negativity
        ##    constraint).  We currently revert to using solve.QP() from
        ##    package quadprog, as constrOptim() already failed in very
        ##    simple test cases.
        ##
        ## Now consider \sum_{p=1}^{k_b} x_{pq}(b).  If k <= k_b for all
        ## b, all M classes from 1 to k are matched to one of the k_b
        ## classes in U(b), hence the sum and also \alpha_q are one.
        ## But then
        ##    \sum_q \bar{m}_{iq}
        ##      =  \sum_b \sum_{p=1}^{k_b} w_b u_{ip}(b) x_{pq}(b)
        ##      <= \sum_b \sum_{p=1}^{k_b} w_b u_{ip}(b)
        ##      =  1
        ## with equality if k = k_b for all b.  I.e.,
        ## * If k = \min_b k_b = \max k_b, we are in case A.
        ## * If k <= \min_b k_b, we are in case B2.
        ## And it makes sense to handle these cases explicitly for
        ## efficiency reasons.

        ## And now for something completely different ... the code.

        k <- .n_of_nonzero_columns(M)
        nr_M <- nrow(M)
        nc_M <- ncol(M)
        nc_memberships <- sapply(memberships, ncol)

        if(k <= min(nc_memberships)) {
            ## Compute the weighted means \bar{M}.
            M <- .weighted_sum_of_matrices(mapply(function(u, p)
                                                  u[ , p[seq_len(k)]],
                                                  memberships,
                                                  permutations,
                                                  SIMPLIFY = FALSE),
                                           w, nr_M)
            ## And add dummy classes if necessary.
            if(k < nc_M)
                M <- cbind(M, matrix(0, nr_M, nc_M - k))
            ## If we always got the same number of classes, we are
            ## done.  Otherwise, equidistribute ...
            if(k < max(nc_memberships))
                M <- pmax(M + (1 - rowSums(M)) / nc_M, 0)
            return(M)
        }

        ## Here comes the general case.

        ## First, compute the \alpha and \beta.
        alpha <- rowSums(rep(w, each = k) *
                         mapply(function(p, n) p[seq_len(k)] <= n,
                                permutations, nc_memberships))
        ## Alternatively (more literally):
        ##   X <- lapply(permutations, .make_X_from_p)
        ##   alpha1 <- double(length = k)
        ##   for(b in seq_along(permutations)) {
        ##       alpha1 <- alpha1 +
        ##           w[b] * colSums(X[[b]][seq_len(nc_memberships[b]), ])
        ##   }
        
        ## A helper function giving suitably permuted memberships.
        pmem <- function(u, p) {
            ## Only matched classes, similar to the one used in value(),
            ## maybe merge eventually ...
            v <- matrix(0, nr_M, k)
            ind <- seq_len(k)
            ind <- ind[p[ind] <= ncol(u)]
            if(any(ind))
                v[ , ind] <- u[ , p[ind]]
            v
        }
        beta <- .weighted_sum_of_matrices(mapply(pmem,
                                                 memberships,
                                                 permutations,
                                                 SIMPLIFY = FALSE),
                                          w, nr_M)
        ## Alternatively (more literally):
        ##   beta1 <- matrix(0, nr_M, nc_M)
        ##   for(b in seq_along(permutations)) {
        ##     ind <- seq_len(nc_memberships[b])
        ##     beta1 <- beta1 +
        ##       w[b] * memberships[[b]][, ind] %*% X[[b]][ind, ]
        ##   }
        
        ## Compute the weighted means \bar{M}.
        M <- .cscale(beta, ifelse(alpha > 0, 1 / alpha, 0))
        ## Alternatively (see comments for .cscale()):
        ##   M1 <- beta %*% diag(ifelse(alpha > 0, 1 / alpha, 0))
        
        ## And add dummy classes if necessary.
        if(k < nc_M)
            M <- cbind(M, matrix(0, nr_M, nc_M - k))

        S <- rowSums(M)
        ## Take care of those rows with row sums < 1.
        ind <- (S < 1)
        if(any(ind)) {
            i_0 <- alpha == 0
            if(any(i_0))
                M[ind, i_0] <- 1 / sum(i_0)
            else
                M[ind, ] <- pmax(M[ind, ] + (1 - S[ind]) / nc_M, 0)
        }
        ## Take care of those rows with row sums > 1.
        ind <- (S > 1)
        if(any(ind)) {
            ## Argh.  Call solve.QP() for each such i.  Alternatively,
            ## could set up on very large QP, but is this any better?
            Dmat <- diag(alpha, nc_M)
            Amat <- t(rbind(rep(-1, nc_M), diag(1, nc_M)))
            bvec <- c(-1, rep(0, nc_M))
            for(i in which(ind))
                M[i, ] <- quadprog::solve.QP(Dmat, alpha * M[i, ],
                                             Amat, bvec)$solution
        }

        M
    }

    memberships <- lapply(clusterings, cl_membership)

    V_opt <- Inf
    M_opt <- NULL
    for(run in seq_along(start)) {
        if(verbose && (nruns > 1L))
            message(gettextf("AOG run: %d", run))
        M <- start[[run]]
        permutations <- lapply(memberships, fit_P, M)
        old_value <- value(M, permutations, memberships, w)
        message(gettextf("Iteration: 0 *** value: %g", old_value))
        iter <- 1L
        while(iter <= maxiter) {
            ## Fit M.
            M <- fit_M(permutations, memberships, w)
            ## Fit \{ P_b \}.
            permutations <- lapply(memberships, fit_P, M)
            ## Update value.
            new_value <- value(M, permutations, memberships, w)
            if(verbose)
                message(gettextf("Iteration: %d *** value: %g",
                                 iter, new_value))
            if(abs(old_value - new_value)
               < reltol * (abs(old_value) + reltol))
                break
            old_value <- new_value
            iter <- iter + 1L
        }
        if(new_value < V_opt) {
            converged <- (iter <= maxiter)
            V_opt <- new_value
            M_opt <- M
        }
        if(verbose)
            message(gettextf("Minimum: %g", V_opt))
    }

    M <- .stochastify(M_opt)
    ## Seems that M is always kept a k columns ... if not, use
    ##   M <- .stochastify(M_opt[, seq_len(k), drop = FALSE])
    rownames(M) <- rownames(memberships[[1L]])
    ## Recompute the value, just making sure ...
    permutations <- lapply(memberships, fit_P, M)
    meta <- list(objval = value(M, permutations, memberships, w),
                 converged = converged)
    M <- .cl_membership_from_memberships(M, k, meta)
    
    as.cl_partition(M)
}

### ** .cl_consensus_partition_GV1

.cl_consensus_partition_GV1 <-
function(clusterings, weights, control)
    .cl_consensus_partition_AOG(clusterings, weights, control, "GV1")


### * .cl_consensus_partition_GV3

.cl_consensus_partition_GV3 <-
function(clusterings, weights, control)
{
    ## Use a SUMT to solve
    ##   \| Y - M M' \|_F^2 => min
    ## where M is a membership matrix and Y = \sum_b w_b M_b M_b'.

    n <- n_of_objects(clusterings)
    max_n_of_classes <- max(sapply(clusterings, n_of_classes))

    ## Control parameters:
    ## k,
    k <- control$k
    if(is.null(k))
        k <- max_n_of_classes
    ## nruns,
    nruns <- control$nruns
    ## start.
    start <- control$start

    w <- weights / sum(weights)

    comemberships <-
        lapply(clusterings, function(x) {
            ## No need to force a common k here.
            tcrossprod(cl_membership(x))
        })

    Y <- .weighted_sum_of_matrices(comemberships, w, n)

    ## Handle start values and number of runs.
    if(!is.null(start)) {
        if(!is.list(start)) {
            ## Be nice to users.
            start <- list(start)
        }
    } else {
        if(is.null(nruns)) {
            ## Use nruns only if start is not given.
            nruns <- 1L
        }
        e <- eigen(Y, symmetric = TRUE)
        ## Use M <- U_k \lambda_k^{1/2}, or random perturbations
        ## thereof.
        M <- e$vectors[, seq_len(k), drop = FALSE] *
            rep(sqrt(e$values[seq_len(k)]), each = n)
        m <- c(M)
        start <- c(list(m),
                   replicate(nruns - 1L,
                             m + rnorm(length(m), sd = sd(m) / sqrt(3)),
                             simplify = FALSE))
    }

    y <- c(Y)

    L <- function(m) sum((y - tcrossprod(matrix(m, n))) ^ 2)
    P <- .make_penalty_function_membership(n, k)
    grad_L <- function(m) {
        M <- matrix(m, n)
        4 * c((tcrossprod(M) - Y) %*% M)
    }
    grad_P <- .make_penalty_gradient_membership(n, k)

    out <- sumt(start, L, P, grad_L, grad_P,
                method = control$method, eps = control$eps,
                q = control$q, verbose = control$verbose,
                control = as.list(control$control))

    M <- .stochastify(matrix(out$x, n))
    rownames(M) <- rownames(cl_membership(clusterings[[1L]]))
    meta <- list(objval = L(c(M)))
    M <- .cl_membership_from_memberships(M, k, meta)

    as.cl_partition(M)
}

### * .cl_consensus_partition_soft_symdiff

.cl_consensus_partition_soft_symdiff <-
function(clusterings, weights, control)
{
    ## Use a SUMT to solve
    ##   \sum_b w_b \sum_{ij} | c_{ij}(b) - c_{ij} | => min
    ## where C(b) = comembership(M(b)) and C = comembership(M) and M is
    ## a membership matrix.

    ## Control parameters:
    ## gradient,
    gradient <- control$gradient
    if(is.null(gradient))
        gradient <- TRUE
    ## k,
    k <- control$k
    ## nruns,
    nruns <- control$nruns
    ## start.
    start <- control$start

    ## Handle start values and number of runs.
    if(!is.null(start)) {
        if(!is.list(start)) {
            ## Be nice to users.
            start <- list(start)
        }
    } else if(is.null(nruns)) {
        ## Use nruns only if start is not given.
        nruns <- 1L
    }
    
    max_n_of_classes <- max(sapply(clusterings, n_of_classes))
    if(is.null(k))
        k <- max_n_of_classes

    B <- length(clusterings)    
    n <- n_of_objects(clusterings)
    
    w <- weights / sum(weights)

    comemberships <-
        lapply(clusterings, function(x) {
            ## No need to force a common k here.
            tcrossprod(cl_membership(x))
        })

    ## Handle start values and number of runs.
    if(!is.null(start)) {
        if(!is.list(start)) {
            ## Be nice to users.
            start <- list(start)
        }
    } else {
        if(is.null(nruns)) {
            ## Use nruns only if start is not given.
            nruns <- 1L
        }
        ## Try using a rank k "root" of the weighted median of the
        ## comemberships as starting value.
        Y <- apply(array(unlist(comemberships), c(n, n, B)), c(1, 2),
                   weighted_median, w)
        e <- eigen(Y, symmetric = TRUE)
        ## Use M <- U_k \lambda_k^{1/2}, or random perturbations
        ## thereof.
        M <- e$vectors[, seq_len(k), drop = FALSE] *
            rep(sqrt(e$values[seq_len(k)]), each = n)
        m <- c(M)
        start <- c(list(m),
                   replicate(nruns - 1L,
                             m + rnorm(length(m), sd = sd(m) / sqrt(3)),
                             simplify = FALSE))
    }

    L <- function(m) {
        M <- matrix(m, n)
        C_M <- tcrossprod(M)
        ## Note that here (as opposed to hard/symdiff) we take soft
        ## partitions as is without replacing them by their closest hard
        ## partitions.
        sum(w * sapply(comemberships,
                       function(C) sum(abs(C_M - C))))
    }
    P <- .make_penalty_function_membership(n, k)
    if(gradient) {
        grad_L <- function(m) {
            M <- matrix(m, n)
            C_M <- tcrossprod(M)
            .weighted_sum_of_matrices(lapply(comemberships,
                                             function(C)
                                             2 * sign(C_M - C) %*% M),
                                      w, n)
        }
        grad_P <- .make_penalty_gradient_membership(n, k)
    }
    else
        grad_L <- grad_P <- NULL

    out <- sumt(start, L, P, grad_L, grad_P,
                method = control$method, eps = control$eps,
                q = control$q, verbose = control$verbose,
                control = as.list(control$control))

    M <- .stochastify(matrix(out$x, n))
    rownames(M) <- rownames(cl_membership(clusterings[[1L]]))
    meta <- list(objval = L(c(M)))
    M <- .cl_membership_from_memberships(M, k, meta)

    as.cl_partition(M)
}

### * .cl_consensus_partition_hard_symdiff

.cl_consensus_partition_hard_symdiff <-
function(clusterings, weights, control)
{
    ## <NOTE>
    ## This is mostly duplicated from relations.
    ## Once this is on CRAN, we could consider having clue suggest
    ## relations ...
    ## </NOTE>
    
    comemberships <-
        lapply(clusterings,
               function(x) {
                   ## Here, we always turn possibly soft partitions to
                   ## their closest hard partitions.
                   ids <- cl_class_ids(x)
                   outer(ids, ids, "==")
                   ## (Simpler than using tcrossprod() on
                   ## cl_membership().)
               })

    ## Could also create a relation ensemble from the comemberships and
    ## call relation_consensus().

    M <- relations:::.make_fit_relation_symdiff_M(comemberships,
                                                  weights)
    k <- control$k
    control <- control$control
    ## Note that currently we provide no support for finding *all*
    ## consensus partitions (but allow for specifying the solver).
    control$all <- FALSE
    I <- if(!is.null(k)) {
        ## <NOTE>
        ## We could actually get the memberships directly in this case.
        relations:::fit_relation_symdiff_E_k(M, k, control)
        ## </NOTE>
    }
    else
        relations:::fit_relation_symdiff(M, "E", control)

    ids <- relations:::get_class_ids_from_incidence(I)
    names(ids) <- cl_object_names(clusterings)

    as.cl_hard_partition(ids)
}
    
### * .cl_consensus_hierarchy_cophenetic

.cl_consensus_hierarchy_cophenetic <-
function(clusterings, weights, control)
{
    ## d <- .weighted_mean_of_object_dissimilarities(clusterings, weights)
    ## Alternatively:
    ## as.cl_dendrogram(ls_fit_ultrametric(d, control = control))
    control <- c(list(weights = weights), control)
    as.cl_dendrogram(ls_fit_ultrametric(clusterings, control = control))
}

### * .cl_consensus_hierarchy_manhattan

.cl_consensus_hierarchy_manhattan <-
function(clusterings, weights, control)
{
    ## Control parameters:
    ## gradient,
    gradient <- control$gradient
    if(is.null(gradient))
        gradient <- TRUE
    ## nruns,
    nruns <- control$nruns
    ## start.
    start <- control$start

    ## Handle start values and number of runs.
    if(!is.null(start)) {
        if(!is.list(start)) {
            ## Be nice to users.
            start <- list(start)
        }
    } else if(is.null(nruns)) {
        ## Use nruns only if start is not given.
        nruns <- 1L
    }
    
    w <- weights / sum(weights)
    B <- length(clusterings)
    ultrametrics <- lapply(clusterings, cl_ultrametric)
    
    if(B == 1L)
        return(as.cl_dendrogram(ultrametrics[[1L]]))

    n <- n_of_objects(ultrametrics[[1L]])
    labels <- cl_object_names(ultrametrics[[1L]])

    ## We need to do
    ##
    ##   \sum_b w_b \sum_{i,j} | u_{ij}(b) - u_{ij} | => min
    ##
    ## over all ultrametrics u.  Let's use a SUMT (for which "gradients"
    ## can optionally be switched off) ...
    
    L <- function(d) {
        sum(w * sapply(ultrametrics, function(u) sum(abs(u - d))))
        ## Could also do something like
        ##   sum(w * sapply(ultrametrics, cl_dissimilarity, d,
        ##                  "manhattan"))
    }
    P <- .make_penalty_function_ultrametric(n)
    if(gradient) {
        grad_L <- function(d) {
            ## "Gradient" is \sum_b w_b sign(d - u(b)).
            .weighted_sum_of_vectors(lapply(ultrametrics,
                                            function(u) sign(d - u)),
                                     w)
        }
        grad_P <- .make_penalty_gradient_ultrametric(n)
    } else
        grad_L <- grad_P <- NULL
            
    if(is.null(start)) {
        ## Initialize by "random shaking" of the weighted median of the
        ## ultrametrics.  Any better ideas?
        ## <NOTE>
        ## Using var(x) / 3 is really L2 ...
        ## </NOTE>
        x <- apply(matrix(unlist(ultrametrics), ncol = B), 1, 
                   weighted_median, w)
        start <- replicate(nruns,  
                           x + rnorm(length(x), sd = sd(x) / sqrt(3)),
                           simplify = FALSE)
    }

    out <- sumt(start, L, P, grad_L, grad_P,
                method = control$method, eps = control$eps,
                q = control$q, verbose = control$verbose,
                control = as.list(control$control))

    d <- .ultrametrify(out$x)
    meta <- list(objval = L(d))
    d <- .cl_ultrametric_from_veclh(d, n, labels, meta)
    
    as.cl_dendrogram(d)
}

### * .cl_consensus_hierarchy_majority

.cl_consensus_hierarchy_majority <-
function(clusterings, weights, control)
{
    w <- weights / sum(weights)

    p <- control$p
    if(is.null(p))
        p <- 1 / 2
    else if(!is.numeric(p) || (length(p) != 1) || (p < 1 / 2) || (p > 1))
        stop("Parameter 'p' must be in [1/2, 1].")

    classes <- lapply(clusterings, cl_classes)
    all_classes <- unique(unlist(classes, recursive = FALSE))
    gamma <- double(length = length(all_classes))
    for(i in seq_along(classes))
        gamma <- gamma + w[i] * !is.na(match(all_classes, classes[[i]]))
    ## Rescale to [0, 1].
    gamma <- gamma / max(gamma)

    maj_classes <- if(p == 1) {
        ## Strict consensus tree.
        all_classes[gamma == 1]
    }
    else
        all_classes[gamma > p]
    attr(maj_classes, "labels") <- attr(classes[[1L]], "labels")

    ## <FIXME>
    ## Stop auto-coercing that to dendrograms once we have suitable ways
    ## of representing n-trees.
    as.cl_hierarchy(.cl_ultrametric_from_classes(maj_classes))
    ## </FIXME>
}

### * Utilities

### ** .cl_consensus_method_default

.cl_consensus_method_default <-
function(type)
{
    switch(type, partition = "SE", hierarchy = "euclidean", NULL)
}

### ** .project_to_leading_columns

.project_to_leading_columns <-
function(x, k)
{
    ## For a given matrix stochastic matrix x, return the stochastic
    ## matrix y which has columns from k+1 on all zero which is closest
    ## to x in the Frobenius distance.
    y <- x[, seq_len(k), drop = FALSE]
    y <- cbind(pmax(y + (1 - rowSums(y)) / k, 0),
               matrix(0, nrow(y), ncol(x) - k))
    ## (Use the pmax to ensure that entries remain nonnegative.)
}

### ** .make_X_from_p

.make_X_from_p <-
function(p) {
    ## X matrix corresponding to permutation p as needed for the AO
    ## algorithms.  I.e., x_{ij} = 1 iff j->p(j)=i.
    X <- matrix(0, length(p), length(p))
    i <- seq_along(p)
    X[cbind(p[i], i)] <- 1
    X
}

### ** .n_of_nonzero_columns

## <NOTE>
## Could turn this into n_of_classes.matrix().
.n_of_nonzero_columns <-
function(x)
    sum(colSums(x) > 0)
## </NOTE>

### ** .cscale

## <FIXME>
## Move to utilities eventually ...
.cscale <-
function(A, x)
{
    ## Scale the columns of matrix A by the elements of vector x.
    ## Formally, A %*% diag(x), but faster.
    ## Could also use sweep(A, 2, x, "*")
    rep(x, each = nrow(A)) * A
}
## </FIXME>

## .make_penalty_function_membership

.make_penalty_function_membership <-
function(nr, nc)
    function(m) {
        sum(pmin(m, 0) ^ 2) + sum((rowSums(matrix(m, nr)) - 1) ^ 2)
    }

## .make_penalty_gradient_membership

.make_penalty_gradient_membership <-
function(nr, nc)
    function(m) {
        2 * (pmin(m, 0) + rep.int(rowSums(matrix(m, nr)) - 1, nc))
    }


### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
