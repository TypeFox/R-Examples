#' Solve a conic optimization problem
#'
#' The function \code{ECOS_csolve} is a wrapper around the ecos
#' \code{csolve} C function. Conic constraints are specified using the
#' \eqn{G} and \eqn{h} parameters and can be \code{NULL} and zero
#' length vector respectively indicating an absence of conic
#' constraints.  Similarly, equality constraints are specified via
#' \eqn{A} and \eqn{b} parameters with \code{NULL} and empty vector
#' values representing a lack of such constraints. At most one of the
#' pair \eqn{(G , h)} or \eqn{(A, b)} is allowed to be absent.
#'
#' @importFrom Matrix sparseMatrix
#'
#' @param c the coefficients of the objective function; the length of this determines the number of variables \eqn{n} in the problem.
#' @param G the inequality constraint sparse matrix in compressed column format, e.g. \link[Matrix]{dgCMatrix-class}. Can be \code{NULL}
#' @param h the right hand size of the inequality constraint. Can be empty numeric vector.
#' @param dims is a list of three named elements: \code{dims['l']} an integer specifying the dimension of positive orthant cone, \code{dims['q']} an integer vector specifying dimensions of second-order cones, \code{dims['e']} an integer specifying the number of exponential cones
#' @param A the optional equality constraint sparse matrix in compressed column format, e.g. \link[Matrix]{dgCMatrix-class}. Can be \code{NULL}
#' @param b the right hand side of the equality constraint, must be specified if \eqn{A} is. Can be empty numeric vector.
#' @param bool_vars the indices of the variables, 1 through \eqn{n}, that are boolean; that is, they are either present or absent in the solution
#' @param int_vars the indices of the variables, 1 through \eqn{n}, that are integers
#' @param control is a named list that controls various optimization parameters; see \link[ECOSolveR]{ecos.control}.
#'
#' @return a list of 8 named items
#'  \describe{
#'   \item{x}{primal variables}
#'   \item{y}{dual variables for equality constraints}
#'   \item{s}{slacks for \eqn{Gx + s <= h}, \eqn{s \in K}}
#'   \item{z}{dual variables for inequality constraints \eqn{s \in K}}
#'   \item{infostring}{gives information about the status of solution}
#'   \item{retcodes}{a named integer vector containing four elements
#'     \describe{
#'       \item{exitflag}{0=\code{OPTIMAL}, 1=\code{PRIMAL INFEASIBLE}, 2=\code{DUAL INFEASIBLE}, -1=\code{MAXIT REACHED}}
#'       \item{iter}{the number of iteration used}
#'       \item{mi_iter}{the number of iterations for mixed integer problems}
#'       \item{numerr}{a non-zero number if a numeric error occurred}
#'     }
#'   }
#'   \item{summary}{a named numeric vector containing
#'     \describe{
#'       \item{pcost}{value of primal objective}
#'       \item{dcost}{value of dual objective}
#'       \item{pres}{primal residual on inequalities and equalities}
#'       \item{dres}{dual residual}
#'       \item{pinf}{primal infeasibility measure}
#'       \item{dinf}{dual infeasibility measure}
#'       \item{pinfres}{primal infeasibility residual}
#'       \item{dinfres}{dual infeasibility residual}
#'       \item{gap}{duality gap}
#'       \item{relgap}{relative duality gap}
#'       \item{r0}{Unknown at the moment to this R package maintainer.}
#'     }
#'   }
#'   \item{timing}{a named numeric vector of timing information consisting of
#'     \describe{
#'       \item{runtime}{the total runtime in ecos}
#'       \item{tsetup}{the time for setup of the problem}
#'       \item{tsolve}{the time to solve the problem}
#'     }
#'   }
#' }
#'
#' @section Details:
#'
#' A call to this function will solve the problem:
#' minimize \eqn{c^Tx}, subject to \eqn{Ax = b}, and \eqn{h - G*x \in K}.
#'
#' Variables can be constrained to be boolean (1 or 0) or integers. This is indicated
#' by specifying parameters \code{bool_vars} and/or \code{int_vars} respectively. If so
#' indicated, the solutions will be found using a branch and bound algorithm.
#'
#' @examples
#'
#' ## githubIssue98
#' G <- local({
#'      Gpr <- c(0.416757847405471, 2.136196095668454, 1.793435585194863, -1.,
#'          0.056266827226329, -1.640270808404989, 0.841747365656204, -1.,
#'          0.416757847405471, 2.136196095668454, 1.793435585194863, -1.,
#'          0.056266827226329, -1.640270808404989, 0.841747365656204, -1., -1.)
#'      Gjc <- as.integer(c(0, 4, 8, 12, 16, 17))
#'      Gir <- as.integer(c(0, 1, 2, 7, 0, 1, 2, 8, 3, 4, 5, 9, 3, 4, 5, 10, 6))
#'      Matrix::sparseMatrix(i = Gir, p = Gjc, x = Gpr, index1 = FALSE)
#' })
#' print(G)
#' c <- as.numeric(c(0, 0, 0, 0, 1))
#' h <- as.numeric(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
#' dims <- list(l = 6L, q = 5L, e = 0L)
#' ECOS_csolve(c = c, G = G, h = h,
#'            dims = dims,
#'            A = NULL, b = numeric(0))
#'
#'
#' ## A larger problem using saved data for the large matrices
#' MPC01 <- readRDS(system.file("misc", "MPC01.rds", package="ECOSolveR"))
#' retval <- ECOS_csolve(c = MPC01$c, G = MPC01$G, h = MPC01$h,
#'                       dims = MPC01$dims)
#' retval$retcodes
#' retval$infostring
#' retval$summary
#'
#' @export
ECOS_csolve <- function(c = numeric(0), G = NULL, h=numeric(0),
                         dims=list(l = integer(0), q = NULL, e = integer(0)),
                         A = NULL, b = numeric(0),
                         bool_vars = integer(0), int_vars = integer(0),
                         control = ecos.control() ) {

    if (!is.null(optionCheck <- checkOptions(control))) {
        stop(optionCheck)
    }

    if (!isNontrivialNumericVector(c)) {
        stop("c should be a nontrivial numeric vector")
    }

    nullG <- is.null(G)
    nontrivialH <- isNontrivialNumericVector(h)

    if ((nullG && nontrivialH) ||
        (!nullG && !nontrivialH)) {
        stop("G and h must be supplied together")
    }

    nullA <- is.null(A)
    nontrivialB <- isNontrivialNumericVector(b)
    if ((nullA && nontrivialB) ||
        (!nullA && !nontrivialB)) {
        stop("A and b must be supplied together")
    }

    if (nullG) {
        m <- 0
        n1 <- length(c)
        Gpr <- Gir <- Gjc <- NULL
    } else {
        if (!inherits(G, "CsparseMatrix")) {
            stop("G is required to be of class dgCMatrix")
        }
        m <- nrow(G)
        n1 <- ncol(G)
        if (m != length(h)) {
            stop("h has incompatible dimension with G")
        }
        Gpr <- G@x
        Gir <- G@i
        Gjc <- G@p
    }

    if (nullA) {
        p <- 0
        n2 <- n1
        Apr <- Air <- Ajc <- NULL
    } else {
        if (!inherits(A, "CsparseMatrix")) {
            stop("A is required to be of class dgCMatrix")
        }
        p <- nrow(A)
        n2 <- ncol(A)
        if (p != length(b)) {
            stop("b has incompatible dimension with A")
        }
        Apr <- A@x
        Air <- A@i
        Ajc <- A@p
    }

    if (n1 != n2) {
        stop("Columns of A and G don't match")
    }


    ## Need to check dims as well
    if (is.null(dims)) {
        stop("dims must be a non-null list")
    }
    ## dimensions of the positive orthant cone
    l <- dims$l
    if (is.null(l)) {
        l <- 0L
    } else {
        if (!isNonnegativeInt(l))
            stop("dims['l'] should be a non-negative int")
    }
    ## dimensions of the second order cones
    q <- dims$q
    if (!is.null(q)) {
        if (typeof(q) != "integer" || !all(q > 0))
            stop("dims['q'] should be an integer vector of positive integers")
    }
    ## number of exponential cones
    e <- dims$e
    if (is.null(e)) {
        e <- 0L
    } else {
        if (!isNonnegativeInt(e))
            stop("dims['e'] should be a non-negative int")
    }
    ## I am not performing this check for now...
    ## check that sum(q) + l = m
    ## if ( (sum(q) + l + 3 * e) != m ) {
    ##     stop("Number of rows of G does not match dims$l + sum(dims$q) + dims$e");
    ## }


    if (typeof(bool_vars) != "integer" || ( length(bool_vars) > 0) && any(bool_vars < 1 | bool_vars > n1) ) {
        stop(sprintf("bool_vars must integers between 1 and %d", n1))
    } else {
        bool_vars <- sort.int(bool_vars - 1L)
    }

    if (typeof(int_vars) != "integer" || ( length(int_vars) > 0) && any(int_vars < 1 | int_vars > n1) ) {
        stop(sprintf("int_vars must integers between 1 and %d", n1))
    } else {
        int_vars <- sort.int(int_vars - 1L)
    }

    ## Need to set up options for the ...
    mnp <- as.integer(c(m, n1, p))


    result <- .Call('ecos_csolve',
                    mnp,
                    l, q, e,
                    Gpr, Gjc, Gir,
                    Apr, Ajc, Air,
                    c, h, b,
                    bool_vars, int_vars,
                    control,
                    package = 'ECOSolveR')
    result
}
