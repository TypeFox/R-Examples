##
##  c o n s t r L s q l i n
##


lsqlincon <- function(C, d,                     # min ||C x - d||_2
                      A = NULL,   b = NULL,     # A x   <= b
                      Aeq = NULL, beq = NULL,   # Aeq x == beq
                      lb = NULL,  ub = NULL)    # lb <= x <= ub
{
    stopifnot(is.numeric(C), is.numeric(d))
    if (is.null(A) && !is.null(b) || !is.null(A) && is.null(b))
        stop("If any, both 'A' and 'b' must be NULL.")
    if (is.null(Aeq) && !is.null(beq) || !is.null(Aeq) && is.null(beq))
        stop("If any, both 'Aeq' and 'beq' must be NULL.")

    if (!is.matrix(C)) C <- matrix(C, 1)
    mc  <- nrow(C);   nc  <- ncol(C);  n <- nc
    if (length(d) != mc)
        stop("Dimensions of 'C' and 'd' do not fit.")
    if (is.null(A) && is.null(Aeq) && is.null(lb) && is.null(ub))
        return(qr.solve(C, d))
    
    if (!is.null(A)) {
        if (!is.matrix(A)) A <- matrix(A, 1)
        ma  <- nrow(A);   na  <- ncol(A)
        if (na != n)
            stop("Number of columns of 'A' does not fit with 'C'.")
        # ATTENTION: quadprog requires  A x >= b !
        A <- -A; b <- -b
    }
    if (is.null(Aeq)) {
        meq <- 0
    } else {
        if (!is.matrix(Aeq)) Aeq <- matrix(Aeq, 1)
        meq  <- nrow(Aeq);   neq  <- ncol(Aeq)
        if (neq != n)
            stop("Number of columns of 'Aeq' does not fit with 'C'.")
    }

    if (is.null(lb)) {
        diag_lb <- NULL
    } else {
        if (length(lb) == 1) {
            lb <- rep(lb, n)
        } else if (length(lb) != n) {
            stop("Length of 'lb' and dimensions of C do not fit.")
        }
        diag_lb <- diag(n)
    }
    if (is.null(ub)) {
        diag_ub <- NULL
    } else {
        if (length(ub) == 1) {
            ub <- rep(ub, n)
        } else if (length(ub) != n) {
            stop("Length of 'ub' and dimensions of C do not fit.")
        }
        # ATTENTION: quadprog requires -x >= -ub
        diag_ub <- -diag(n)
        ub <- -ub
    }
    
    Dmat <- t(C) %*% C
    dvec <- t(C) %*% d
    
    Amat <- rbind(Aeq, A, diag_lb, diag_ub)
    bvec <- c(beq, b, lb, ub)

    rslt <- quadprog::solve.QP(Dmat, dvec, t(Amat), bvec, meq=meq)
    rslt$solution
}

