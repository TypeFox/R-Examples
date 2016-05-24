#' Deterministic g/Q functions - examples and templates
#' 
#' Template for the \code{deterministic.g.function} 
#' argument to \code{\link{ltmle}} or \code{\link{ltmleMSM}}.
#' 
#' \code{MaintainTreatment} and \code{MaintainControl} are two commonly used
#' \code{deterministic.g.function}s.
#' 
#' The intended use of the templates is for the user to copy and paste the
#' function arguments and body and then fill in the required sections. They
#' will not run as-is. Note that there are no comments in the functions as
#' saved. Versions with comments may be found in Examples section below.
#' 
#' MaintainTreatment and MaintainControl may be passed as-is for the
#' \code{deterministic.g.function} argument to \code{\link{ltmle}} or
#' \code{\link{ltmleMSM}}
#' 
#' Note that censoring nodes in \code{data} may be passed as binaries but they
#' are converted to the preferred format of factors with levels "censored" and
#' "uncensored" before deterministic functions are called.  Also note that
#' nodes may be passed to ltmle as either the names of nodes or numerical
#' column indicies, but they are all converted to numerical indicies before
#' deterministic functions are called.  If the \code{survivalFunction} argument
#' to \code{ltmle} or \code{ltmleMSM} is \code{TRUE}, the package automatically
#' assumes that once Y jumps to 1, all future Y nodes stay 1 and treatment does
#' not change. It is not necessary to specify this in deterministic functions.
#' 
#' @aliases MaintainTreatment MaintainControl deterministic.g.function_template
#' deterministic.Q.function_template
#' @param data the 'data' data.frame passed to \code{ltmle} or \code{ltmleMSM}
#' @param current.node the column index of data corresponding to the A or C
#' node (for g) or L or Y node (for Q)
#' @param nodes list of column indicies, components: \itemize{ \item \code{A}
#' Anodes (treatment) \item \code{C} Cnodes (censoring) \item \code{L} Lnodes
#' (time-varying covariates) \item \code{Y} Ynodes (events) \item \code{AC}
#' Anodes and Cnodes combined and sorted \item \code{LY} Lnodes and Ynodes
#' combined, sorted, "blocks" removed - see \code{\link{ltmle}} }
#' @param called.from.estimate.g TRUE or FALSE - your function will be called
#' with \code{called.from.estimate.g=TRUE} during estimation of g and
#' \code{called.from.estimate.g=FALSE} during estimation of Q.
#' @return A deterministic.g.function should return a list with components:
#' \item{is.deterministic }{vector of logicals, length=nrow(data)} \item{prob1
#' }{the probability that data[is.deterministic, current.node] == 1, vector of
#' length 1 or length(which(is.deterministic))} A deterministic.Q.function
#' should return a list with components: \item{is.deterministic }{vector of
#' logicals, length=nrow(data)} \item{Q.value}{the iterated expectation of the
#' final Y, vector of length 1 or length(which(is.deterministic))}
#' 
#' NOTE: The \code{Q.value} component is not used or required when
#' \code{called.from.estimate.g} is \code{TRUE}
#' @author Joshua Schwab \email{joshuaschwab@@yahoo.com}
#' @seealso \code{\link{ltmle}}, \code{\link{ltmleMSM}}
#' @examples
#' 
#' # Show template for a deterministic.g.function (comments will not be
#' # shown, see below for comments)
#' deterministic.g.function_template
#' 
#' # Show template for a deterministic.Q.function (comments will not be
#' # shown, see below for comments)
#' deterministic.Q.function_template
#' 
#' # Use MaintainTreatment
#' set.seed(1)
#' rexpit <- function(x) rbinom(n = length(x), size = 1, prob = plogis(x))
#' n <- 100
#' W <- rnorm(n)
#' A1 <- rexpit(W)
#' A2 <- as.numeric(rexpit(W) | A1)  #treatment at time 1 implies treatment at time 2
#' Y <- rexpit(W + A1 + A2 + rnorm(n))
#' data <- data.frame(W, A1, A2, Y)
#' 
#' result <- ltmle(data, Anodes = c("A1", "A2"), Ynodes = "Y", abar = c(1, 1), 
#'     deterministic.g.function = MaintainTreatment)
#' 
#' # deterministic.g.function_template with comments:
#' 
#' deterministic.g.function_template <- function(data, current.node, nodes) {
#'     # data: the 'data' data.frame passed to ltmle/ltmleMSM current.node: the
#'     # column index of data corresponding to the A or C node (see
#'     # is.deterministic below) nodes: list of column indicies, components: A,
#'     # C, L, Y, AC (Anodes and Cnodes combined and sorted), LY (Lnodes and
#'     # Ynodes combined, sorted, 'blocks' removed - see ?ltmle) Note that nodes
#'     # may be passed to ltmle as either the names of nodes or numerical column
#'     # indicies, but they are all converted to numerical indicies before
#'     # deterministic.g.function is called
#'     
#'     # deterministic.g.function will be called at all Anodes and Cnodes
#'     # return(NULL) is equivalent to return(list(is.deterministic=rep(FALSE,
#'     # nrow(data)), prob1=numeric(0)))
#'     
#'     # define is.deterministic here: vector of logicals, length=nrow(data)
#'     # define prob1 here: the probability that data[is.deterministic,
#'     # current.node] == 1, vector of length 1 or
#'     # length(which(is.deterministic))
#'     is.deterministic <- stop("replace me!")
#'     prob1 <- stop("replace me!")
#'     return(list(is.deterministic = is.deterministic, prob1 = prob1))
#' }
#' 
#' # deterministic.Q.function_template with comments:
#' 
#' deterministic.Q.function_template <- function(data, current.node, nodes, 
#'     called.from.estimate.g) {
#'     # data: the 'data' data.frame passed to ltmle/ltmleMSM current.node: the
#'     # column index of data corresponding to the A or C node (see
#'     # is.deterministic below) nodes: list of column indicies, components: A,
#'     # C, L, Y, AC (Anodes and Cnodes combined and sorted), LY (Lnodes and
#'     # Ynodes combined, sorted, 'blocks' removed - see ?ltmle)
#'     # called.from.estimate.g: TRUE or FALSE - your function will be called
#'     # with called.from.estimate.g=TRUE during estimation of g and
#'     # called.from.estimate.g=FALSE during estimation of Q. During estimation
#'     # of g, only the is.deterministic element of the return list will be
#'     # used.  Note that nodes may be passed to ltmle as either the names of
#'     # nodes or numerical column indicies, but they are all converted to
#'     # numerical indicies before deterministic.Q.function is called
#'     
#'     # It is not necessary to specify that deterministic Y events (Y==1)
#'     # indicate a deterministic Q value of 1; this is automatic.
#'     # deterministic.Q.function will be called at all Lnodes and Ynodes (after
#'     # removing 'blocks') and Anodes and Cnodes (see called.from.estimate.g
#'     # above) return(NULL) is equivalent to
#'     # return(list(is.deterministic=rep(FALSE, nrow(data)),
#'     # Q.value=numeric(0)))
#'     
#'     # define is.deterministic here: vector of logicals, length=nrow(data)
#'     # define Q.value here: the iterated expectation of the final Y, vector of
#'     # length 1 or length(which(is.deterministic))
#'     is.deterministic <- stop("replace me!")
#'     Q.value <- stop("replace me!")
#'     return(list(is.deterministic = is.deterministic, Q.value = Q.value))
#' }
#' 
#' @export deterministic.g.function_template
deterministic.g.function_template <- function(data, current.node, nodes) {
  # data: the 'data' data.frame passed to ltmle/ltmleMSM
  # current.node: the column index of data corresponding to the A or C node (see is.deterministic below)
  # nodes: list of column indicies, components: A, C, L, Y, AC (Anodes and Cnodes combined and sorted), 
  #   LY (Lnodes and Ynodes combined, sorted, "blocks" removed - see ?ltmle)
  # Note that nodes may be passed to ltmle as either the names of nodes or numerical column indicies, but they
  #   are all converted to numerical indicies before deterministic.g.function is called
  
  # deterministic.g.function will be called at all Anodes and Cnodes
  # return(NULL) is equivalent to return(list(is.deterministic=rep(FALSE, nrow(data)), prob1=numeric(0)))
  
  #define is.deterministic here: vector of logicals, length=nrow(data)
  #define prob1 here: the probability that data[is.deterministic, current.node] == 1, 
  #  vector of length 1 or length(which(is.deterministic))
  is.deterministic <- stop("replace me!")
  prob1 <- stop("replace me!")
  return(list(is.deterministic=is.deterministic, prob1=prob1))  
}

#' @describeIn deterministic.g.function_template Template for the \code{deterministic.Q.function} 
#' argument to \code{\link{ltmle}} or \code{\link{ltmleMSM}}.
#' @export
deterministic.Q.function_template <- function(data, current.node, nodes, called.from.estimate.g) {
  # data: the 'data' data.frame passed to ltmle/ltmleMSM
  # current.node: the column index of data corresponding to the A or C node (see is.deterministic below)
  # nodes: list of column indicies, components: A, C, L, Y, AC (Anodes and Cnodes combined and sorted), 
  #   LY (Lnodes and Ynodes combined, sorted, "blocks" removed - see ?ltmle)
  # called.from.estimate.g: TRUE or FALSE - your function will be called with called.from.estimate.g=TRUE during 
  #   estimation of g and called.from.estimate.g=FALSE during estimation of Q. During estimation of g, only
  #   the is.deterministic element of the return list will be used.
  # Note that nodes may be passed to ltmle as either the names of nodes or numerical column indicies, but they
  #   are all converted to numerical indicies before deterministic.Q.function is called
  
  # It is not necessary to specify that deterministic Y events (Y==1) indicate 
  #   a deterministic Q value of 1; this is automatic.
  # deterministic.Q.function will be called at all Lnodes and Ynodes (after removing "blocks") 
  #   and Anodes and Cnodes (see called.from.estimate.g above)
  # return(NULL) is equivalent to return(list(is.deterministic=rep(FALSE, nrow(data)), Q.value=numeric(0)))
  
  #define is.deterministic here: vector of logicals, length=nrow(data)
  #define Q.value here: the iterated expectation of the final Y, vector of length 1 or length(which(is.deterministic))
  is.deterministic <- stop("replace me!")
  Q.value <- stop("replace me!")
  return(list(is.deterministic=is.deterministic, Q.value=Q.value))  
}

#' @export
MaintainTreatment <- function(data, current.node, nodes) {
  #if the previous Anode is 1, all subsequent Anodes are 1 
  Anodes <- nodes$A
  if (!(current.node %in% Anodes)) return(NULL)
  if (!(any(Anodes < current.node))) return(NULL)
  
  prev.a.node <- max(Anodes[Anodes < current.node])
  is.deterministic <- data[, prev.a.node] == 1
  return(list(is.deterministic=is.deterministic, prob1=1))  
}

#' @export
MaintainControl <- function(data, current.node, nodes) {
  #if the previous Anode is 0, all subsequent Anodes are 0
  Anodes <- nodes$A
  if (!(current.node %in% Anodes)) return(NULL)
  if (!(any(Anodes < current.node))) return(NULL)
  
  prev.a.node <- max(Anodes[Anodes < current.node])
  is.deterministic <- data[, prev.a.node] == 0
  return(list(is.deterministic=is.deterministic, prob1=0))  
}
