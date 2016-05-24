#' Manipulate the component parts of formulas, expressions, calls, name/symbols 
#' and list and vectors of such objects.
#' 
#' \code{lhs, rhs, op, and op.type } retrieve the various parts of R formulas, 
#' calls, expressions, names/symbols.  These functions were designed to 
#' greatly facilitate symbolic manupulation using native R objects.  Also 
#' provided are methods to handle list of these objects. 
#' 
#' @param x object from where to get/set the lhs/rhs
#' @param value the value to set for the lhs/rhs
#' @param ... arguments passed to additional methods
#' 
#' \code{lhs} retrieves the left-hand side 
#' \code{rhs} retrieves the right-hand side 
#' \code{op}  retrieves the operation 
#' \code{op.type} returns the type operator
#'
#' There are also functions \code{lhs.vars} and \code{rhs.vars}.  Like 
#' \code{ \link{all.vars} }, these functions interpret the variables on the 
#' left-hand and right-hand sides respectively.
#' 
#' These are simple functions for extracting the left-hand side, right-hand 
#' side, operator and operator type from formulas, expressions, calls, 
#' names/symbols and list containing these objects.  lhs, rhs are only defined 
#' for formulas and calls ( and list and expressions ) that are defined with 
#' either one of the relational or tilde ('~') operators. If the object does 
#' not contain one of these operators, it will fail with a warning.
#'
#' The defined operator types are defined by the operator.tools package: See 
#' \code{\link[operator.tools]{operators}} and
#' \code{\link[operator.tools]{setOperator}}
#'
#' The \code{lhs.vars} and \code{rhs.vars} methods, return the variables used on
#' the lhs and rhs, respectively.  If special formula variables are used, such 
#' as '.', a data.frame or environment must also be provided such that the 
#' variable list may be properly infered.
#' 
#' @return Value depends on the argument.
#' 
#' @author Christopher Brown
#' 
#' @seealso terms, all.vars, all.names, \code{\link[operator.tools]{operators}}
#' 
#' @examples 
#' 
#'   # FORMULA
#'   f <- A + B ~ C + D
#'   lhs(f)
#'   lhs(f) <- quote( E / F )
#'
#'   rhs(f)
#'   rhs(f) <- quote( G + H ) 
#'   op(f)
#'   op(rhs(f))
#'   op( quote(A) )  # NULL: 
#'   op.type(f)
#'
#'   # ONE-SIDED FORMULA
#'   f <- ~ A   # 
#'   lhs(f)     # NULL
#'   rhs(f)     # A
#'
#'
#'   # EXPRESSION
#'   e <- expression( A + B == C + D )
#'   lhs(e)
#'   rhs(e)
#'   op(e)
#'   op.type(e)
#'
#'
#'   # CALL
#'   c <- quote( A + B > C + D )
#'   lhs(c)
#'   lhs(c) <- quote(E)
#'   rhs(c)
#'
#'   op(c)
#'   op.type(c)
#'
#'   # ASSIGNMENT 
#'   a  <- quote( A <- B ) 
#'   lhs(a)
#'   rhs(a) 
#'   op(a)
#'   op.type(a) 
#'
#' @name formula.parts
#' @rdname formula.parts
#' @docType methods
#' @import operator.tools
#' @import methods

NULL
