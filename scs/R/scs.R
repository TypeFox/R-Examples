#  ---------------------------------------------------------
#  scs
#  ===
#' @title SCS - Splitting Conic Solver 
#'
#' @description Solves convex cone programs via operator splitting. 
#' @param A a matrix of constraint coefficients. 
#'        \bold{NOTE:} The rows of matrix A have to be ordered according to the 
#'        order given in subsection \dQuote{Allowed cone parameters}. For more information see \bold{README}.
#' @param b a numeric vector giving the primal constraints
#' @param obj a numeric vector giving the primal objective
#' @param cone a list giving the cone sizes
#' @param control a list giving the control parameters. For more information see \bold{README}.
#' @return list of solution vectors x, y, s and information about run
#' @details 
#' 
#' A more detailed description can be found in the README file.
#' 
#' \subsection{Important Note}{ \cr
#' The order of the rows in matrix \eqn{A} has to correspond to the order given in 
#' the table \dQuote{Cone Arguments}, which means means rows corresponding to 
#' \emph{primal zero cones} should be first, rows corresponding to \emph{non-negative cones} second, 
#' rows corresponding to \emph{second-order cone} third, rows corresponding to \emph{positive semidefinite cones} fourth, 
#' rows corresponding to \emph{exponential cones} fifth and rows corresponding to \emph{power cones} at last.
#' }
#'
#' \subsection{SCS can solve}{
#' \enumerate{
#'   \item linear programs (LPs)
#'   \item second-order cone programs (SOCPs)
#'   \item semidefinite programs (SDPs)
#'   \item exponential cone programs (ECPs)
#'   \item power cone programs (PCPs)
#'   \item problems with any combination of cones, which can be defined by the parameters listed in the subsection \dQuote{Allowed cone parameters}
#' } }
#'
#' \subsection{Allowed \emph{cone} parameters are}{
#' \tabular{rllll}{ 
#'    \tab \bold{Parameter} \tab \bold{Type} \tab \bold{Length} \tab             \bold{Description}                       \cr
#'    \tab \code{f}         \tab integer     \tab \eqn{1}       \tab number of primal zero cones (dual free cones),       \cr
#'    \tab                  \tab             \tab               \tab which corresponds to the primal equality constraints \cr
#'    \tab \code{l}         \tab integer     \tab \eqn{1}       \tab number of linear cones (non-negative cones)          \cr
#'    \tab \code{q}         \tab integer     \tab \eqn{\geq1}   \tab vector of second-order cone sizes                    \cr
#'    \tab \code{s}         \tab integer     \tab \eqn{\geq1}   \tab vector of positive semidefinite cone sizes           \cr
#'    \tab \code{ep}        \tab integer     \tab \eqn{1}       \tab number of primal exponential cones                   \cr
#'    \tab \code{ed}        \tab integer     \tab \eqn{1}       \tab number of dual exponential cones                     \cr
#'    \tab \code{p}         \tab numeric     \tab \eqn{\geq1}   \tab vector of primal/dual power cone parameters          
#' } }
#'
#' \subsection{Allowed \emph{control} parameters are}{
#' \tabular{rllll}{ 
#'    \tab \bold{Parameter} \tab \bold{Type}    \tab             \bold{Description}                                          \tab \bold{Default} \cr
#'    \tab \code{max_iters} \tab integer        \tab giving the maximum number of iterations                                 \tab   2500   \cr
#'    \tab \code{normalize} \tab boolean        \tab heuristic data rescaling                                                \tab   TRUE   \cr
#'    \tab \code{verbose}   \tab boolean        \tab write out progress                                                      \tab   FALSE  \cr
#'    \tab \code{cg_rate}   \tab numeric        \tab for indirect, tolerance goes down like \eqn{\frac{1}{iter}^{cg\_rate}}  \tab      2   \cr
#'    \tab \code{scale}     \tab numeric        \tab if normalized, rescales by this factor                                  \tab      5   \cr
#'    \tab \code{rho_x}     \tab numeric        \tab x equality constraint scaling                                           \tab   1e-3   \cr
#'    \tab \code{alpha}     \tab numeric        \tab relaxation parameter                                                    \tab    1.5   \cr
#'    \tab \code{eps}       \tab numeric        \tab convergence tolerance                                                   \tab   1e-3
#' } }
#' @examples
#' A <- matrix(c(1, 1), ncol=1)
#' b <- c(1, 1)
#' obj <- 1
#' cone <- list(f = 2)
#' control <- list(eps = 1e-3, max_iters = 50)
#' sol <- scs(A, b, obj, cone, control) 
#' sol
#  max_iters=2500L, normalize=TRUE, verbose=TRUE,
#  cg_rate=2.0, scale=5.0, rho_x=1e-03, alpha=1.5, eps=1e-3
#  ---------------------------------------------------------
scs <- function(A, b, obj, cone, control=list(max_iters=2500L, normalize=TRUE, verbose=FALSE,
                cg_rate=2.0, scale=5.0, rho_x=1e-03, alpha=1.5, eps=1e-6)) {

    if ( class(A) == "simple_triplet_matrix" ) {
        ## sparseMatrix returns an object of class "dgCMatrix"
        A <- sparseMatrix( i=A$i, j=A$j, x=A$v, dims=c(A$nrow, A$ncol) )
    } else {
        A <- as(A, "dgCMatrix")
    }

    data <- list(m = dim(A)[1], n = dim(A)[2], 
                 Ax = A@x, Ai = A@i, Ap = A@p, b = b, c = obj)
    
    ret <- .Call("scsr", data, cone, control, PACKAGE="scs")
    return(ret)
}

