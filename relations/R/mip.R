## A simple framework for representing and solving mixed integer linear
## programs (MILPs) of the form
##   optimize obj' x
##   such that mat %*% x dir rhs
## and mixed integer quadratic programs (MIQPs) of the form
##   optimize x' Q x / 2 + c' x
##   such that mat %*% x dir rhs
## with possibly given types (C/I/B for continuous/integer/binary) and
## given additional (lower and upper) bounds on x.
## (Default of course x >= 0).

### * MILPs

MILP <-
function(objective, constraints, bounds = NULL, types = NULL,
         maximum = FALSE)
{
    ## Currently, 'constraints' always is a (not necessarily named) list
    ## with mat, dir and rhs, which really is the most general case of
    ## linear constraints we can think of.  Let us add names for now;
    ## eventually, there should be more sanity checking and maybe a
    ## creator for linear constraint objects.

    names(constraints) <- c("mat", "dir", "rhs")
    .structure(list(objective = objective, constraints = constraints,
                    bounds = bounds, types = types, maximum = maximum),
               class = "MILP")
}

.MILP_solvers <-
    c("glpk", "lpsolve", "symphony", "cplex")

solve_MILP <-
function(x, solver = NULL, control = list())
{
    ## <NOTE>
    ## Ideally, we would use some registration mechanism for solvers.
    ## Currently, there is only little support for control arguments.
    ## In particular, one cannot directly pass arguments to the solver.
    ## </NOTE>

    ## Handle the boundary case of no variables.
    if(!length(x$objective)) {
        y <- .solve_empty_MIP(x)
        if(!is.null(nos <- control$n)
           && !identical(as.integer(nos), 1L))
            y <- list(y)
        return(y)
    }

    solver <- match.arg(solver, .MILP_solvers)

    ## If more than one (binary) solution is sought and the solver does
    ## not provide direct support, use poor person's branch and cut:
    if(!is.null(nos <- control$n) && (solver != "cplex")) {
        control$n <- NULL
        ## Mimic the mechanism currently employed by Rcplex(): return a
        ## list of solutions only if nos > 1 (or NA).
        if(!identical(as.integer(nos), 1L)) {
            add <- identical(control$add, TRUE)
            control$add <- NULL
            return(.find_up_to_n_binary_MILP_solutions(x, nos, add,
                                                       solver, control))
        }
    }
    ## Note that lpSolve could find all binary solutions for all-binary
    ## programs.

    switch(solver,
           "cplex" = .solve_MILP_via_cplex(x, control),
           "glpk" = .solve_MILP_via_glpk(x),
           "lpsolve" = .solve_MILP_via_lpsolve(x),
           "symphony" = .solve_MILP_via_symphony(x))
}

### * MIQPs

MIQP <-
function(objective, constraints, bounds = NULL, types = NULL,
         maximum = FALSE)
{
    ## See MILP() for the comments on constraint objects.
    ## We also add names to the list of objective coefficients.
    names(objective) <- c("Q", "L")
    names(constraints) <- c("mat", "dir", "rhs")
    .structure(list(objective = objective, constraints = constraints,
                    bounds = bounds, types = types, maximum = maximum),
               class = "MIQP")
}

.MIQP_solvers <-
    c("cplex")

solve_MIQP <-
function(x, solver = NULL, control = list())
{
    ## Currently, only CPLEX can generally be used for solving MIQPs.
    ## For the other MILP solvers, all-binary programs can be solved via
    ## linearization.
    ## <NOTE>
    ## Actually, linearization only requires that the quadratic part is
    ## all-binary.  Maybe support the mixed linear part case eventually.
    ## </NOTE>

    ## Handle the boundary case of no variables.
    if(!length(x$objective)) {
        y <- .solve_empty_MIP(x)
        if(!is.null(nos <- control$n)
           && !identical(as.integer(nos), 1L))
            y <- list(y)
        return(y)
    }

    is_BQP <- identical(unique(x$types), "B")
    solver <- if(is_BQP) {
        ## If this is an all-binary problem, we can linearize, so use
        ## non-commercial defaults (obviously, there should eventually
        ## be a way to specify the default MILP and MIQP solver).
        match.arg(solver, .MILP_solvers)
    } else {
        ## Use a MIQP solver by default.
        match.arg(solver, c(.MIQP_solvers, .MILP_solvers))
    }
    ## For real MIQP solvers (currently only CPLEX), do not linearize by
    ## default, but allow for doing so for debugging purposes.
    if(solver %in% .MIQP_solvers) {
        if(identical(control$linearize, TRUE))
            .solve_BQP_via_linearization(x, solver, control)
        else {
            ## Add switch() when adding support for other MIQP solvers.
            .solve_MIQP_via_cplex(x, control)
        }
    } else {
        ## If this is an all-binary problem, we can linearize.
        if(is_BQP)
            .solve_BQP_via_linearization(x, solver, control)
        else
            stop(gettextf("Solver '%s' can only handle all-binary quadratic programs.",
                          solver),
                 domain = NA)
    }
}

.solve_BQP_via_linearization <-
function(x, solver, control)
{
    ## Number of variables.
    n <- length(x$objective$L)
    ## Solve via linearization.
    y <- solve_MILP(.linearize_BQP(x), solver, control)
    ## Reduce solution to the original variables.
    finisher <- function(e) {
        e$solution <- e$solution[seq_len(n)]
        e
    }
    ## <FIXME>
    ## Wouldn't it be simpler to check if y inherits from MIP_solution?
    if(!is.null(nos <- control$n) && !identical(nos, 1L))
        lapply(y, finisher)
    else
        finisher(y)
    ## </FIXME>
}

### * Solver interfaces

### ** CPLEX

.solve_MILP_via_cplex <-
function(x, control)
{
    ## Wrap into the common MIP CPLEX framework.
    x$objective <- list(Q = NULL, L = x$objective)
    .solve_MIP_via_cplex(x, control)
}

.solve_MIQP_via_cplex <-
function(x, control)
{
    ## Ensure that the coefficient matrix of the quadratic term is
    ## symmetric, as required by Rcplex.
    Q <- x$objective$Q
    ## <CHECK>
    ## Does Rcplex really support simple triplet matrix Q coefficients?
    ## If not, would need
    ##   x$objective$Q <- as.matrix((Q + t(Q)) / 2)
    ## instead:
    x$objective$Q <- (Q + t(Q)) / 2
    ## </CHECK>
    .solve_MIP_via_cplex(x, control)
}

.solve_MIP_via_cplex <-
function(x, control)
{
    ## Currently, no direct support for bounds.
    ## <FIXME>
    ## Should expand the given bounds and map into lb/ub arguments.
    if(!is.null(x$bounds))
        stop("Solver currently does not support variable bounds.")
    ## </FIXME>

    .as_Rcplex_sense <- function(x) {
        TABLE <- c("L", "L", "G", "G", "E")
        names(TABLE) <- c("<", "<=", ">", ">=", "==")
        TABLE[x]
    }
    sense <- .as_Rcplex_sense(x$constraints$dir)

    types <- .expand_types(x$types, length(x$objective$L))

    mat <- x$constraints$mat
    if(is.simple_triplet_matrix(mat)) {
        ## Reorder indices as CPLEX needs a column major order
        ## representation i.e., column indices j have to be in ascending
        ## order.
        column_major_order <- order(mat$j)
        mat$i <- mat$i[column_major_order]
        mat$j <- mat$j[column_major_order]
        mat$v <- mat$v[column_major_order]
    } else {
        mat <- as.matrix(mat)
    }

    if(is.null(nos <- control$n)) nos <- 1L
    value_is_list_of_solutions <- !identical(as.integer(nos), 1L)

    out <-
        tryCatch(Rcplex::Rcplex(Qmat = x$objective$Q,
                                cvec = x$objective$L,
                                Amat = mat,
                                sense = sense,
                                bvec = x$constraints$rhs,
                                vtype = types,
                                objsense = if(x$maximum) "max" else "min",
                                control = list(trace = 0, round = 1),
                                n = nos
                                ),
                 error = identity)
    if(inherits(out, "error")) {
        ## Explicitly catch and rethrow CPLEX unavailability errors.
        msg <- conditionMessage(out)
        if(regexpr("Could not open CPLEX environment\\.", msg) > -1L)
           stop(msg, call. = FALSE)
        ## Currently, Rcplex signals problems via error() rather than
        ## returning a non-zero status.  Hence, we try catching these
        ## errors.  (Of course, these could also be real errors ...).
        solution <- rep(NA_real_, length(types))
        objval <- NA_real_
        status <- 2                     # or whatever ...
        names(status) <- msg            # should be of length one ...
        out <- .make_MIP_solution(solution, objval, status)
        if(value_is_list_of_solutions) out <- list(out)
    } else {
        out <- if(value_is_list_of_solutions)
            lapply(out, .canonicalize_solution_from_cplex, x)
        else
            .canonicalize_solution_from_cplex(out, x)
    }
    out
}

.canonicalize_solution_from_cplex <-
function(out, x)
{
    solution <- out$xopt
    ## For the time being ...
    ## Since Rcplex 0.1-4 integers are rounded (via control argument
    ## 'round' which we set accordingly when calling Rcplex()) but no
    ## new optimal solution based on these values is calculated.  Hence,
    ## we no longer round ourselves, but recompute objval.
    objval <- sum(solution * x$objective$L)
    if(!is.null(Q <- x$objective$Q))
        objval <- objval + .xtQx(Q, solution) / 2
    status <- out$status
    ## Simple db for "ok" status results:
    ok_status_db <-
        c("CPX_STAT_OPTIMAL" = 1L,      # (Simplex or barrier): optimal
                                        # solution is available
          "CPXMIP_OPTIMAL" = 101L,      # (MIP): optimal integer solution
                                        # has been found
          "CPXMIP_OPTIMAL_TOL" = 102L,  # (MIP): Optimal soluton with
                                        # the tolerance defined by epgap
                                        # or epagap has been found
          "CPXMIP_POPULATESOL_LIM" = 128L, # (MIP-MultSols): The limit on
                                        # mixed integer solutions
                                        # generated by populate has been
                                        # reached
          "CPXMIP_OPTIMAL_POPULATED" = 129L, # (MIP-MultSols): Populate
                                        # has completed the enumeration of
                                        # all solutions it could enumerate
          "CPXMIP_OPTIMAL_POPULATED_TOL" = 130L # (MIP-MultSols): similar
                                        # to 129L but additionally
                                        # objective value fits the
                                        # tolerance specified by paramaters
          )
    status <- ifelse(status %in% ok_status_db, 0, status)
    .make_MIP_solution(solution, objval, status)
}

### ** GLPK

.solve_MILP_via_glpk <-
function(x)
{
    out <- Rglpk::Rglpk_solve_LP(x$objective,
                                 x$constraints$mat,
                                 x$constraints$dir,
                                 x$constraints$rhs,
                                 bounds = x$bounds,
                                 types = x$types,
                                 max = x$maximum)
    .make_MIP_solution(out$solution, out$optimum, out$status)
}

### ** lp_solve

.solve_MILP_via_lpsolve <-
function(x)
{
    ## Currently, no direct support for bounds.
    ## <FIXME>
    ## Should rewrite the given bounds into additional constraints.
    if(!is.null(x$bounds))
        stop("Solver currently does not support variable bounds.")
    ## </FIXME>

    types <- .expand_types(x$types, length(x$objective))

    ## Version 5.6.1 of lpSolve has added sparse matrix support via
    ## formal 'dense.const' as well as binary variable types.
    mat <- x$constraints$mat
    out <- if(is.simple_triplet_matrix(mat)) {
        ## In the sparse case, lpSolve currently (2008-11-22) checks
        ## that every constraint is used in the sense that each row has
        ## at least one entry.  So if for some reason this is not the
        ## case, let us add one zero entry for such rows (note that we
        ## cannot simply drop them as the corresponding constraint may
        ## be violated).
        ind <- which(tabulate(mat$i, mat$nrow) == 0)
        if(len <- length(ind)) {
            mat$i <- c(mat$i, ind)
            mat$j <- c(mat$j, rep.int(1L, len))
            mat$v <- c(mat$v, rep.int(0, len))
        }
        lpSolve::lp(if(x$maximum) "max" else "min",
                    x$objective,
                    const.dir = x$constraints$dir,
                    const.rhs = x$constraints$rhs,
                    int.vec = which(types == "I"),
                    binary.vec = which(types == "B"),
                    dense.const = cbind(mat$i, mat$j, mat$v))
    } else {
        lpSolve::lp(if(x$maximum) "max" else "min",
                    x$objective,
                    as.matrix(mat),
                    x$constraints$dir,
                    x$constraints$rhs,
                    int.vec = which(types == "I"),
                    binary.vec = which(types == "B"))
    }
    status_db <-
        ## Solver status values from lp_lib.h:
       c("UNKNOWNERROR" = -5L,
         "DATAIGNORED" = -4L,
         "NOBFP" = -3L,
         "NOMEMORY" = -2L,
         "NOTRUN" = -1L,
         "OPTIMAL" = 0L,
         "SUBOPTIMAL" = 1L,
         "INFEASIBLE" = 2L,
         "UNBOUNDED" = 3L,
         "DEGENERATE" = 4L,
         "NUMFAILURE" = 5L,
         "USERABORT" = 6L,
         "TIMEOUT" = 7L,
         "RUNNING" = 8L,
         "PRESOLVED" = 9L)
    status <- status_db[match(out$status, status_db)]
    solution <- out$solution
    objval <- if(status == 0L)
        sum(solution * out$objective)
    else
        out$objval
    .make_MIP_solution(solution, objval, status)
}

### ** SYMPHONY

.solve_MILP_via_symphony <-
function(x)
{
    out <- Rsymphony::Rsymphony_solve_LP(x$objective,
                                         x$constraints$mat,
                                         x$constraints$dir,
                                         x$constraints$rhs,
                                         bounds = x$bounds,
                                         types = x$types,
                                         max = x$maximum)
    .make_MIP_solution(out$solution, out$objval, out$status)
}

### * Utilities

.expand_types <-
function(x, n)
{
    if(is.null(x)) {
        ## Continuous by default.
        rep.int("C", n)
    }
    else {
        if(!is.character(x) || !all(x %in% c("C", "I", "B")))
            stop("Invalid MIP variable types.")
        ## Be nicer than necessary ...
        rep(x, length.out = n)
    }
}

.find_up_to_n_binary_MILP_solutions <-
function(x, nos = 1L, add = FALSE, solver = NULL, control = NULL)
{
    ## Find up to n binary MILP solutions using a simple branch and cut
    ## approach (repeatedly splitting the binary variables and cutting
    ## the non-optimal branches).

    if(is.na(nos))
        nos <- .Machine$integer.max

    y <- solve_MILP(x, solver, control)
    if((y$status != 0) || (nos == 1L) || !any(x$types == "B"))
        return(list(y))

    Vopt <- y$objval
    tol <- 1e-8
    ## (Smaller than .Machine$double.eps^0.5 as e.g. used in the
    ## all.equal() comparisons.)
    ## We used to have 1e-10, but SYMPHONY was not as precise as this.
    v_is_not_optimal <- if(x$maximum)
        function(v) v < Vopt - tol
    else
        function(v) v > Vopt + tol
    ## Improved versions could use relative tolerance, and/or take the
    ## number of variables into account.  Or maybe we can figure out the
    ## tolerance employed by the solvers when they declare optimality?

    if(add) {
        ## Find up to n solutions by adding binary constraints for each
        ## split.  This is most space efficient, but typically takes
        ## considerably more time that the default based on successive
        ## reductions by substitution.
        return(.find_up_to_n_binary_MILP_solutions_via_add(x, nos,
                                                           solver,
                                                           control,
                                                           y,
                                                           v_is_not_optimal))
    }

    ## Find solutions by recursively splitting binary variable and
    ## substituting the values into the program.

    ## Suppose x_j is one of the binary variables.  Let x[-j] be the
    ## vector of variables after dropping x_j.  If x_j has value b_j,
    ## the objective function is
    ##   c[-j]' x[-j] + c_j b_j,
    ## and the constraints become
    ##   (mat[, -j] %*% x[-j]) dir (rhs - mat[, j] * b_j)
    ## Note that the new constraint matrix and dir do not depend on the
    ## value of b_j.
    ## When recursively splitting, we need to keep track of the b
    ## values, the sum of the c_j b_j terms to add to the reduced
    ## objective value to obtain the value of the original problem, and
    ## the right hand sides.
    ## Also, binary positions are a nuisance to keep track of, so we
    ## start by rearranging variables to have the binary variables come
    ## last in reverse order of splitting.

    n_of_variables <- length(x$objective)
    types <- .expand_types(x$types, n_of_variables)
    binary_positions <- which(types == "B")
    n_of_binary_variables <- length(binary_positions)

    verbose <- identical(control$verbose, TRUE)

    .make_node <- function(b, y, v, r) list(b = b, y = y, v = v, r = r)

    .split_single_binary_variable <-
        function(node, i, pos_i, obj_i, mat_i, pos_b) {
            ## Try to avoid unnecessary copying of x via <<-
            ## manipulations.
            ## We know that the solution with the current b[i] is
            ## optimal, so we try the effect of flipping b[i].
            b <- node$b
            node$y$solution <- node$y$solution[-pos_i]
            if(b[i] == 0) {
                b1 <- b
                b1[i] <- 1
                v1 <- node$v + obj_i
                r1 <- node$r - mat_i
                x$constraints$rhs <<- r1
                y1 <- solve_MILP(x, solver, control)
                ## Uncomment for debugging ...
                ## if(verbose) {
                ##     V <- y1$objval + v1
                ##     split <- (y1$status == 0) && !v_is_not_optimal(V)
                ##     message(sprintf("b[i] = 0, flip objval: %f, Delta: %.12f, status: %d, split: %s",
                ##                     V, Vopt - V, y1$status, split))
                ## }
                if((y1$status != 0) || v_is_not_optimal(y1$objval + v1))
                    list(node)
                else {
                    ## Fill the rest of b with the optimal entries.
                    b1[-seq_len(i)] <- y1$solution[pos_b]
                    list(node, .make_node(b1, y1, v1, r1))
                }
            } else {
                b0 <- b
                b0[i] <- 0
                v0 <- node$v
                node$v <- v0 + obj_i
                r0 <- node$r
                node$r <- r0 - mat_i
                x$constraints$rhs <<- r0
                y0 <- solve_MILP(x, solver, control)
                ## Uncomment for debugging ...
                ## if(verbose) {
                ##     V <- y0$objval + v0
                ##     split <- (y0$status == 0) && !v_is_not_optimal(V)
                ##     message(sprintf("b[i] = 1, flip objval: %f, Delta: %.12f, status: %d, split: %s",
                ##                     V, Vopt - V, y0$status, split))
                ## }
                if((y0$status != 0) || v_is_not_optimal(y0$objval + v0))
                    list(node)
                else {
                    ## Fill the rest of b with the optimal entries.
                    b0[-seq_len(i)] <- y0$solution[pos_b]
                    list(node, .make_node(b0, y0, v0, r0))
                }
            }
        }

    ## We allow callers to specify the order in which binary splits
    ## should be attempted (so the i-th split is for binary variable
    ## order[i]).
    order <- control$order
    if(is.null(order) || length(order) != n_of_binary_variables)
        order <- seq_len(n_of_binary_variables)
    ## Rearrange variables to have the binary ones last in reverse split
    ## order.
    ind <- c(which(types != "B"), rev(which(types == "B")[order]))
    x$objective <- x$objective[ind]
    x$constraints$mat <- x$constraints$mat[, ind, drop = FALSE]
    x$types <- x$types[ind]
    y$solution <- y$solution[ind]

    pos_i <- n_of_variables
    pos_b <- seq(from = n_of_variables, length.out =
                 n_of_binary_variables, by = -1L)
    nodes <- list(.make_node(y$solution[pos_b], y, 0, x$constraints$rhs))
    for(i in seq_len(n_of_binary_variables)) {
        pos_b <- pos_b[-1L]
        obj_i <- x$objective[pos_i]
        mat_i <- c(x$constraints$mat[, pos_i])
        x$objective <- x$objective[-pos_i]
        x$constraints$mat <- x$constraints$mat[, -pos_i, drop = FALSE]
        x$types <- x$types[-pos_i]
        nodes <- do.call(c,
                         lapply(nodes,
                                .split_single_binary_variable,
                                i, pos_i, obj_i, mat_i, pos_b))
        len <- length(nodes)
        if(verbose)
            message(gettextf("N_of_binary_variables: %d *** N_of_optimal_branches: %d",
                             i, len))
        if(len >= nos) {
            nodes <- nodes[seq_len(nos)]
            break
        }
        pos_i <- pos_i - 1L
    }

    pos <- order(c(which(types != "B"), which(types == "B")[order]))
    finisher <- function(node) {
        ## Need to reconstruct solutions from the binary and non-binary
        ## parts and the respective objective values.
        y <- node$y
        y$solution <- c(y$solution, node$b)[pos]
        y$objval <- y$objval + node$v
        y
    }

    lapply(nodes, finisher)
}

.find_up_to_n_binary_MILP_solutions_via_add <-
function(x, nos = 1L, solver = NULL, control = NULL, y,
         v_is_not_optimal)
{
    ## Recursively add 0/1 constraints for the binary variables.
    ## Note that one can do this via adding the additional constraints
    ## to the original ones, or by maintaining them separately (e.g., in
    ## a vector of values for the binary variables).  We do the latter
    ## which uses as little space as possible, but requires extra time
    ## for merging both constraints when solving augmented problems.
    ## Alternatively, one could allow choosing between either approach.
    ## Note also that we need to keep new constraints and corresponding
    ## solutions to avoid recomputing solutions when enough were found.

    n_of_variables <- length(x$objective)
    types <- .expand_types(x$types, n_of_variables)
    binary_positions <- which(types == "B")
    n_of_binary_variables <- length(binary_positions)

    verbose <- identical(control$verbose, TRUE)

    mat <- x$constraints$mat

    .make_additional_constraint_matrix <-
        if(is.simple_triplet_matrix(mat))
            function(bpos) {
                len <- length(bpos)
                simple_triplet_matrix(seq_len(len), bpos,
                                      rep.int(1, len),
                                      len, n_of_variables)
            }
        else
            function(bpos) {
                len <- length(bpos)
                add <- matrix(0, len, n_of_variables)
                add[cbind(seq_len(len), bpos)] <- 1
                add
            }

    .solve_MILP_with_additional_binary_constraints <-
        function(x, bpos, bval) {
            len <- length(bpos)
            x$constraints <-
                list(mat = rbind(mat,
                                 .make_additional_constraint_matrix(bpos)),
                     dir = c(x$constraints$dir, rep.int("==", len)),
                     rhs = c(x$constraints$rhs, bval))
            solve_MILP(x, solver, control)
        }

    .split_single_binary_variable <- function(y, i) {
        oi <- order[i]
        ind <- order[seq_len(i)]
        pos <- binary_positions[ind]
        b <- y$solution[binary_positions]
        ## Try flipping the i-th binary variable and see whether this
        ## also delivers optimal solutions.
        b[oi] <- 1 - b[oi]
        yf <- .solve_MILP_with_additional_binary_constraints(x, pos, b[ind])
        if((yf$status != 0) || v_is_not_optimal(yf$objval))
            list(y)
        else
            list(y, yf)
    }

    ylist <- list(y)
    ## We allow callers to specify the order in which binary splits
    ## should be attempted (so the i-th split is for binary variable
    ## order[i], i.e., at binary_positions[order[i]]).
    order <- control$order
    if(is.null(order))
        order <- seq_len(n_of_binary_variables)
    for(i in seq_len(n_of_binary_variables)) {
        ylist <-
            do.call(c, lapply(ylist, .split_single_binary_variable, i))
        len <- length(ylist)
        if(verbose)
            message(gettextf("N_of_binary_variables: %d *** N_of_optimal_branches: %d",
                             i, len))
        if(len >= nos) {
            ylist <- ylist[seq_len(nos)]
            break
        }
    }

    ylist
}

.linearize_BQP <-
function(x)
{
    ## Linearize an all-binary quadratic program
    ##   \sum_{i,j} q_{ij} x_i x_j / 2 + \sum_i c_i x_i
    ## as described e.g. in "Pseudo-Boolean Optimization" by E. Boros
    ## and P. Hammer (boros01pseudoboolean.pdf): rewrite the criterion
    ## function as
    ##   \sum_{i < j} r_{ij} y_{ij} + \sum_i s_i x_i
    ## with
    ##  r_{ij} = (q_{ij} + q_{ji}) / 2
    ##  s_i = c_i + q_{ii} / 2
    ## and the additional constraints
    ##   y_{ij} <= x_i, y_{ij} <= x_j              (A)
    ##   y_{ij} >= 0, y_{ij} >= x_i + x_j - 1      (B)
    ## where for a minimization problem (A) is redundant if r_{ij} > 0
    ## and (B) if \gamma_{ij} < 0, and vice versa for a maximization
    ## problem.

    if(!inherits(x, "MIQP") && !identical(unique(x$types), "B"))
        stop("Can only linearize all-binary quadratic programs.")

    ## Could do some sanity checking here.
    Q <- x$objective$Q
    c <- x$objective$L
    n <- length(c)
    R <- (Q + t(Q)) / 2
    if(is.simple_triplet_matrix(Q)) {
        ## Transform coefficients.
        ## Cannot easily have a diag() method for simple triplet
        ## matrices.
        s <- c + Q[cbind(seq_len(n), seq_len(n))] / 2
        ## Quadratic coefficients and respective variables.
        p <- (R$i < R$j) & (R$v != 0)
        i <- R$i[p]
        j <- R$j[p]
        r <- R$v[p]
    } else {
        ## Transform coefficients.
        s <- c + diag(Q) / 2
        ## Quadratic coefficients and respective variables.
        I <- upper.tri(R)
        r <- R[I]
        p <- which(r != 0)
        I <- which(I, arr.ind = TRUE)
        i <- I[p, 1L]
        j <- I[p, 2L]
        r <- r[p]
    }
    nr <- length(r)

    ## Constraints.
    mat <- x$constraints$mat
    pn <- which(r < 0)                  # Negative positions.
    pp <- which(r > 0)                  # Positive positions.
    ## <NOTE>
    ## To experiment with not dropping redundant constraints, do:
    ##    pn <- pp <- seq_along(r)
    ## </NOTE>
    npn <- length(pn)
    npp <- length(pp)
    if(x$maximum) {
        if(is.simple_triplet_matrix(mat)) {
            add_i <- c(rep.int(seq_len(npp), 2L),
                       rep.int(seq_len(npp) + npp, 2L),
                       rep.int(seq_len(npn) + 2L * npp, 3L))
            add_j <- c(i[pp], n + pp,
                       j[pp], n + pp,
                       i[pn], j[pn], n + pn)
            add_v <- rep.int(c(-1, 1, -1, 1, -1, 1),
                             c(npp, npp, npp, npp, 2L * npn, npn))
            mat <- rbind(cbind(mat,
                               simple_triplet_zero_matrix(nrow(mat), nr)),
                         simple_triplet_matrix(add_i, add_j, add_v,
                                               npn + 2L * npp, n + nr))
        } else {
            add <- matrix(0, npn + 2L * npp, n + nr)
            ## Constraints
            ##    y_{ij} <= x_i, y_{ij} <= x_j             (A)
            ## if r_{ij} > 0:
            ind <- seq_len(npp)
            add[cbind(ind, i[pp])] <- -1
            add[cbind(ind, n + pp)] <- 1
            ind <- ind + npp
            add[cbind(ind, j[pp])] <- -1
            add[cbind(ind, n + pp)] <- 1
            ## Constraints
            ##   y_{ij} >= 0, y_{ij} >= x_i + x_j - 1      (B)
            ## if r_{ij} < 0 (where the former is implicit):
            ind <- seq_len(npn) + 2L * npp
            add[cbind(ind, i[pn])] <- -1
            add[cbind(ind, j[pn])] <- -1
            add[cbind(ind, n + pn)] <- 1
            mat <- rbind(cbind(mat, matrix(0, nrow(mat), nr)), add)
        }
        dir <- c(x$constraints$dir,
                 rep.int("<=", 2L * npp),
                 rep.int(">=", npn))
        rhs <- c(x$constraints$rhs,
                 rep.int(0, 2L * npp),
                 rep.int(-1, npn))
    } else {
        if(is.simple_triplet_matrix(mat)) {
            add_i <- c(rep.int(seq_len(npn), 2L),
                       rep.int(seq_len(npn) + npn, 2L),
                       rep.int(seq_len(npp) + 2L * npn, 3L))
            add_j <- c(i[pn], n + pn,
                       j[pn], n + pn,
                       i[pp], j[pp], n + pp)
            add_v <- rep.int(c(-1, 1, -1, 1, -1, 1),
                             c(npn, npn, npn, npn, 2L * npp, npp))
            mat <- rbind(cbind(mat,
                               simple_triplet_zero_matrix(nrow(mat), nr)),
                         simple_triplet_matrix(add_i, add_j, add_v,
                                               npp + 2L * npn, n + nr))
        } else {
            add <- matrix(0, 2L * npn + npp, n + nr)
            ## Constraints
            ##    y_{ij} <= x_i, y_{ij} <= x_j             (A)
            ## if r_{ij} < 0:
            ind <- seq_len(npn)
            add[cbind(ind, i[pn])] <- -1
            add[cbind(ind, n + pn)] <- 1
            ind <- ind + npn
            add[cbind(ind, j[pn])] <- -1
            add[cbind(ind, n + pn)] <- 1
            ## Constraints
            ##   y_{ij} >= 0, y_{ij} >= x_i + x_j - 1      (B)
            ## if r_{ij} > 0 (where the former is implicit):
            ind <- seq_len(npp) + 2L * npn
            add[cbind(ind, i[pp])] <- -1
            add[cbind(ind, j[pp])] <- -1
            add[cbind(ind, n + pp)] <- 1
            mat <- rbind(cbind(mat, matrix(0, nrow(mat), nr)), add)
            dir <- c(x$constraints$dir,
                     rep.int("<=", 2L * npn),
                     rep.int(">=", npp))
            rhs <- c(x$constraints$rhs,
                     rep.int(0, 2L * npn),
                     rep.int(-1, npp))
        }
    }

    MILP(c(s, r),
         list(mat, dir, rhs),
         x$bounds,
         rep.int(c("B", "C"), c(n, nr)),
         x$maximum)
}

.make_MIP_solution <-
function(solution, objval, status, ...)
    .structure(list(solution = solution,
                    objval = objval,
                    status = status, ...),
               class = "MIP_solution")

.make_types <-
function(n, I = NULL, B = NULL)
{
    ## Create MIP variable types spec from possibly given positions of
    ## integer and binary variables.
    types <- rep.int("C", n)
    if(!is.null(I)) types[I] <- "I"
    if(!is.null(B)) types[B] <- "B"
    types
}

.relax_mixed_integer_program <-
function(x)
{
    ## Relax MILP or MIQP by dropping integrality constraints (I -> C),
    ## and changing binary constraints x \in \{0, 1\} to x \in [0,1].
    if(!inherits(x, "MILP") && !inherits(x, "MIQP"))
        stop("Can only relax mixed integer linear or quadratic programs.")

    mat <- x$constraints$mat
    n_of_variables <- ncol(mat)
    types <- .expand_types(x$types, n_of_variables)
    binary_positions <- which(types == "B")
    if(n_of_binary_variables <- length(binary_positions)) {
        ## For binary variables x_i, we need to add the constraint
        ## x_i <= 1.
        if(is.simple_triplet_matrix(mat))
            add <-
                simple_triplet_matrix(seq_len(n_of_binary_variables),
                                      binary_positions,
                                      rep.int(1, n_of_binary_variables),
                                      n_of_binary_variables,
                                      n_of_variables)
        else {
            add <- matrix(0, n_of_binary_variables, n_of_variables)
            add[cbind(seq_len(n_of_binary_variables),
                      binary_positions)] <- 1
        }
        x$constraints$mat <- rbind(mat, add)
        x$constraints$dir <- c(x$constraints$dir,
                               rep.int("<=", n_of_binary_variables))
        x$constraints$rhs <- c(x$constraints$rhs,
                               rep.int(1, n_of_binary_variables))
    }
    x$types <- rep.int("C", n_of_variables)
    x
}

.solve_empty_MIP <-
function(x)
{
    ## Check whether constraints are satisfied (interpreting each lhs as
    ## empty sum with value 0):
    constraints <- split(x$constraints$rhs, x$constraints$dir)
    if(all(unlist(Map(function(dir, rhs) get(dir)(0, rhs),
                      names(constraints), constraints))))
        .make_MIP_solution(double(), 0, 0L)
    else
        .make_MIP_solution(double(), NA_real_, 2L)
}

.xtQx <-
function(Q, x)
{
    ## Value of quadratic form t(x) %*% Q %*% x.
    ## As we implement simple triplet matrices in S3, we could only have
    ## %*% and crossprod methods if we created S3 generics for these ...
    if(is.simple_triplet_matrix(Q))
        sum(Q$v * x[Q$i] * x[Q$j])
    else
        c(crossprod(x, Q %*% x))
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
