##
##  l i n p r o g . R  Linear Programming Solver
##


linprog <- function(cc, A = NULL, b = NULL, Aeq = NULL, beq = NULL,
                    lb = NULL, ub = NULL, x0 = NULL, I0 = NULL,
                    bigM = 100, maxiter = 20, maximize = FALSE)
{
    .lp.check(cc, A, b, Aeq, beq, lb, ub, x0, I0)
    if (maximize) cc <- -cc

    # there are only inequality constraints
    if (is.null(Aeq) && all(b >= 0)) {
        # go and find a base vector and index set
        Pin <- .lp.create(cc, A, b, lb, ub)
        Sin <- .lp.solve(Pin$cc, Pin$A, Pin$b, Pin$x0, Pin$I0,
                         maxiter = maxiter)

    # find a feasible base vector first
    } else {
        # there are only equality constraints
        if (is.null(A) && is.null(lb) && is.null(ub)) {
            A <- Aeq
            b <- beq
            cc0 <- cc

        # there are both equality and inequality constraints
        } else {
            Pin <- .lp.create(cc, A, b, lb, ub)
            Ain <- Pin$A
            bin <- Pin$b
            cc0 <- Pin$cc
            if (!is.null(Aeq)) {
                Aeq <- cbind(Aeq, matrix(0, nrow = nrow(Aeq),
                                    ncol = ncol(Ain)-ncol(Aeq)))
                A <- rbind(Ain, Aeq)
                b <- c(Pin$b, beq)
            } else {
                A <- Ain
                b <- bin
            }
        }

        # only equalities; now guarantee that b >= 0
        inds <- which(b < 0)
        if (length(inds) > 0) {
            A[inds, ] <- -A[inds, ]
            b[inds] <- -b[inds]

        }
        # all constraints combined in one equality matrix A and vectors b, cc
        m <- nrow(A); n <- ncol(A)
        A0 <- cbind(A, diag(m))
        b0 <- b

        x0 <- c(rep(0, n), b0)
        I0 <- (n+1):(n+m)

        if (bigM > 0) {                 # big-M approach
            M  <- bigM
            k <- 1; kmax <- 8
            while (k < kmax) {
                f <- c(cc0, rep(M, m))
                Sin <- .lp.solve(f, A0, b0, x0, I0, maxiter = maxiter)
                if (Sin$errno < 0 || all(Sin$x[I0] == 0)) break
                M <- 10 * M
                k <- k + 1
            }
            # browser()
            if (k >= kmax) {
                Sin$errno <- -4
            }
        } else {                        # phase-I approach
            f <- rep(1, n+m)
            Sin <- .lp.solve(f, A0, b0, x0, I0, maxiter = maxiter)
            if (Sin$errno < 0) {
                Sin$errno <- -5
            } else {
                x1 <- Sin$x[1:n]; z1 <- Sin$x[(n+1):(n+m)]
                I1 <- which(x1 != 0)
                if (length(I1) != m) {
                    Sin$errno <- -5
                } else if (all(z1 == 0)) {
                    Sin <- .lp.solve(cc0, A, b, x1, I1, maxiter = maxiter)
                } else
                    Sin$errno <- -2
            }
        }
    }
    errno <- Sin$errno
    errmsg <- c("Maximum no. of iterations exceeded.",
                "Problem is most likely unfeasible.",
                "Problem is most likely unbounded.",
                "Big-M appraoch not successful.",
                "Phase-I approach not successful.")

    if (errno < 0) {
        mess <- errmsg[abs(errno)]
        return(list(x = NA, fval = NA, errno = errno, message = mess))
    } else {
        mess <- "Solver LP converged successfully."
        x <- Sin$x[1:length(cc)]
        y <- Sin$y[1:length(b)]
        f <- sum(cc * x)
        if (maximize) f <- -f
        return(list(x = x,  fval = f, errno = errno, message = mess))
    }
    cat("LP solver stopped.\n")
}


.lp.check <- function(cc, A, b, Aeq, beq, lb, ub, x0, I0) {
    stopifnot(is.numeric(cc))
    if (is.numeric(lb) && length(lb) != length(cc) ||
        is.numeric(ub) && length(ub) != length(cc))
        stop("Arguments 'cc', 'lb', and 'ub' must have equal lengths.")

    if (!is.null(A)) {
        if (is.null(b)) stop("Vector 'b' must be present if 'A' is.")
        if (!is.matrix(A)) stop("Argument 'A' must be a matrix.")
        m <- nrow(A); n <- ncol(A)
        if (length(cc) != n) stop("Length of 'cc' does not fit to 'ncol(A)'.")
        if (length(b) != m)  stop("Length of 'b' does not fit to 'nrow(A)'.")
    }
}


.lp.create <- function(cc, Ain, bin, lb, ub) {
    A <- Ain
    b <- bin
    if (is.null(Ain)) {
        m_in <- 0
        n <- length(cc)
    } else {
        m_in <- nrow(Ain)
        n <- ncol(Ain)
    }

    # upper bounds
    if (!is.null(ub)) {
        inds <- which(is.finite(ub))
        m1 <- length(inds)
        B <- matrix(0, nrow = m1, ncol = n)
        for (k in 1:m1) B[k, inds[k]] <- 1
        A <- rbind(A, B)
        b <- c(b, ub[inds])
    }

    # lower bounds
    if (!is.null(lb)) {
        inds <- which(lb > 0)   # all(lb >= 0)
        m2 <- length(inds)
        B <- matrix(0, nrow = m2, ncol = n)
        for (k in 1:m2) B[k, inds[k]] <- -1
        A <- rbind(A, B)
        b <- c(b, -lb[inds])
    }

    # slack variables and base vector
    m <- nrow(A)            # m <- m_in + m1 + m2
    A <- cbind(A, diag(m))
    cc <- c(cc, rep(0, m))
    v <- c(rep(0, n), b)
    I <- (n+1):(n+m)

    return(list(A = A, b = b, cc = cc, x0 = v, I0 = I))
}


.lp.solve <- function(cc, A, b, x, I, maxiter = maxiter) {
    # Linear program  A x = b, minimize cc * x.
    # x is assumed to be a base vector with index set I.
    stopifnot(is.numeric(A), is.numeric(b), is.numeric(x))

    if (!is.matrix(A))
        stop("Argument 'A' must be a numeric matrix.")
    m <- nrow(A); n <- ncol(A)
    if (length(b) != m || length(x) != n)
        stop("One of arguments 'x' or 'b' does not have correct length.")

    # cat("Start value:", cc %*% x, "\n")
    # check I: must be a strict subset of {1, ..., n}

    iter <- 0
    errno <- 0
    while (iter < maxiter) {
        iter <- iter + 1
        J <- setdiff(1:n, I)
        B <- A[, I]; N <- A[, J]
        x_I <- x[I]; x_J <- x[J]
        c_I <- cc[I]; c_J <- cc[J]
        
        # Transform base vector
        y <- solve(t(B), c_I)
        u <- numeric(n)
        for (j in J)
	        u[j] <- cc[j] - t(A[, j]) %*% y
        
        # Case 1: x is a solution of the LP
        if (all(u >= 0)) {
            # cat("Solution found!\n")
            errno <- 1
            break
        }
        
        # Case 2: r in J with u[j] < 0
        r <- which(u < 0)
        if (length(r) > 1) r <- r[sample(1:length(r), 1)]
        
        dd <- solve(B, A[, r])
        d <- numeric(n)
        d[I] <- dd
        
        # Case 2a: The LP does not have a solution
        if (all(d <= 0)) {
            # cat("LP is infeasible!\n")
            errno <- -2
            break
        }
        
        # Case 2b: i in I with d[i] > 0
        i0 <- which(d > 0)
        t <- min(x[i0]/d[i0])
        s <- which(x[i0]/d[i0] == t); s <- i0[s]
        if (length(s) > 1)
            s <- s[sample(1:length(s), 1)]
        
        z <- numeric(n)
        z[r] <- t
        z[I] <- x_I - t*d[I]

        # Define a new base vector
        x <- z
        is <- which(I == s)
        I[is] <- r
        # cat("New base vector:", cc %*% x, "\n")
    }
    if (iter >= maxiter) {
        errno <- -1
    }
    if (errno < 0) {
        x <- y <- NA
    }

    return(list(x = x, dual.x = y, errno = errno))
}
