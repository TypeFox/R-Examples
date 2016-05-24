## ROI plugin: quadprog
## based on quadprog package

## BASIC SOLVER METHOD
solve_QP <- function( x, control ) {
    ## quadprog does not support variable bounds per se, thus we add
    ## them as constraints
    x <- ROI:::as.no_V_bounds_OP( x )

    ## since ROI 0.1 the objective function is an STM, quadprog only supports 1 (dense) row
    L <- if( slam::is.simple_triplet_matrix(terms(objective(x))$L) ){
        stopifnot( dim(terms(objective(x))$L)[1] == 1L )
        as.numeric( as.matrix(terms(objective(x))$L) )
    } else {
        terms(objective(x))$L
    }

    ## quadprog needs an appropiately formated linear part of the objective function
    if( !length(L) )
        L <- double(length(objective(x)))
    stopifnot( length(L) == length(objective(x)) )

    ## solve the QP
    out <- .quadprog_solve_QP(Q   = terms(objective(x))$Q,
                              L   = L,
                              mat = constraints(x)$L,
                              dir = constraints(x)$dir,
                              rhs = constraints(x)$rhs,
                              max = x$maximum,
                              control = control)
    ROI:::canonicalize_solution( solution = out$sol,
                                 optimum  = objective(x)(out$sol),
                                 status   = out$ierr,
                                 solver   = ROI:::get_solver_name(getPackageName()),
                                 solver_return_object = out)
}

## SOLVER SUBMETHODS
.quadprog_solve_QP <- function(Q, L, mat, dir, rhs, max, control) {

  ## Partially borrowed from the fPortfolio function '.rquadprog'
  ## Description:
  ##   Goldfarb and Idnani's quadprog solver function
  ## Note:
  ##   Requires to load contributed R package quadprog from which we
  ##   directly use the Fortran subroutine of the quadratic solver.
  ## Note 2: we currently disable the direct call of the Fortran routine
  ##   since R-devel 2013-11-27 requires a different mechanism in calling .Fortran
  ##   (see 5.4 and 5.4.2 in Writing R Extensions)

  ind_eq  <- which( dir == "==")
  ind_geq <- which( (dir == ">=") | (dir == ">") )
  ind_leq <- which( (dir == "<=") | (dir == "<") )
  meq <- length(ind_eq)

  ## everything except equaltity constraints have to be ">="
  ## FIXME: no replacement method for 'simple_triplet_matrix[i, ]<-'
  ##       thus, coercing to 'matrix'
  Amat <- as.matrix(mat)
  Amat[ ind_leq, ] <- -Amat[ ind_leq, ]
  bvec <- rhs
  bvec[ ind_leq ] <- -bvec[ ind_leq ]

  ## We have to sort constraints. The first meq constraints are
  ## equality constraints
  if( length(ind_eq) ) {
    Amat <- rbind( Amat[ ind_eq, ], Amat[ -ind_eq, ] )
    bvec <- c( bvec[ ind_eq ], bvec[ -ind_eq ] )
  }
  ## quadprog uses mat^T in the constraints
  Amat <- t(Amat)
  ## replace Inf with .Machine$double.xmax
  Amat[ is.infinite(Amat) & (Amat <= 0) ] <- -.Machine$double.xmax
  Amat[ is.infinite(Amat) & (Amat >= 0) ] <-  .Machine$double.xmax
  bvec[ is.infinite(bvec) & (bvec <= 0) ] <- -.Machine$double.xmax
  bvec[ is.infinite(bvec) & (bvec >= 0) ] <-  .Machine$double.xmax

  ## dvec in objective function according to direction of optimization
  dvec <- if( max )
    L
  else
    -L

  ## same with the quadratic term
  Dmat <- if( max )
    -as.matrix(Q)
  else
    as.matrix(Q)

  factorized <- control$factorized
  if( is.null(factorized) )
    factorized <- formals(solve.QP)$factorized

  ## FIXME: temporarily added in order to remove .Fortran call
  ##########  ##########  ##########  ##########  ##########
  out <- tryCatch( solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = meq, factorized = factorized), error = identity )
  out$ierr <- 0L
  if( inherits(out, "error") ){
    ierr <- if( out$message == "constraints are inconsistent, no solution!" )
      1L
    else if (out$message == "matrix D in quadratic function is not positive definite!" )
      2L
    else
      3L ## new status code
    out <- list( solution = rep(NA, nrow(Dmat)), ierr = ierr )
  }
  out$sol <- out$solution # we need this since the Fortran routine would provide this field
  out
  ## END FIXME: temporarily added in order to remove .Fortran call
  ##########  ##########  ##########  ##########  ##########

  ## ## number objectives
  ## n_obj <- nrow(Dmat)
  ## ## number constraints
  ## n_constr <- ncol(Amat)

  ## r = min(n_obj, n_constr)
  ## work = rep(0, 2 * n_obj+ r * (r + 5)/2 + 2 * n_constr + 1)

  ## ## FIXME: do we need santiy checks here?

  ## # Optimize:
  ## .Fortran(.QP_qpgen2,
  ##          as.double(Dmat),
  ##          dvec = as.double(dvec),
  ##          as.integer(n_obj),
  ##          as.integer(n_obj),
  ##          sol = double(n_obj),
  ##          lagr = double(n_constr),
  ##          crval = double(1),
  ##          as.double(Amat),
  ##          as.double(bvec),
  ##          as.integer(n_obj),
  ##          as.integer(n_constr),
  ##          as.integer(meq),
  ##          iact = integer(n_constr),
  ##          nact = 0L,
  ##          iter = integer(2L),
  ##          work = as.double(work),
  ##         ierr = 0L) # , NAOK = TRUE

}


## STATUS CODES
.add_status_codes <- function(){
    ## quadprog
    solver <- ROI:::get_solver_name( getPackageName() )
    ROI:::add_status_code_to_db(solver,
                                0L,
                                "OPTIMAL",
                                "Solution is optimal",
                                0L
                                )
    ROI:::add_status_code_to_db(solver,
                                1L,
                                "INCONSISTENT",
                                "Constraints are inconsistent, no solution."
                                )
    ROI:::add_status_code_to_db(solver,
                                2L,
                                "NOT_POSITIVE_DEFINITE",
                                "quadratic term in function is not positive definite."
                                )
    ## FIXME: temporary status code until Fortran routine is called directly again
    ROI:::add_status_code_to_db(solver,
                                3L,
                                "ROI_INTERFACE_ERROR",
                                "contact the plugin maintainer."
                                )
    invisible(TRUE)
}
