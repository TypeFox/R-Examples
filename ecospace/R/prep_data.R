#' Internal Ecospace Functions.
#'
#' Internal functions not intended to be called directly by users.
#'
#' @param ecospace An ecospace framework (functional trait space) of class
#'   \code{ecospace}.
#' @param Smax Maximum number of species (or other taxa) to include in
#'   simulation.
#'
#' @details Pre-allocate the basic data frame structure so classes (and factor
#'   levels, if used) are properly set. Bit clunky but efficient solution:
#'   adding first placeholder column allows \code{cbind} to work properly with
#'   data frame later, but requires deleting that placeholder at end.
#'
#' @return Returns an emptry data frame with pre-defined number of species and
#'   functional traits of proper class type.
#'
#' @seealso \code{\link{create_ecospace}} for how to create an ecospace
#'   framework (functional trait space).
#'
#' @export
prep_data <- function(ecospace, Smax) {
  if(class(ecospace)!="ecospace") stop("use create.ecospace() to create a formal ecospace framework.\n")
  pre.data <- as.data.frame(factor(Smax))
  nchar <- length(ecospace) - 1
  seq <- seq_len(nchar)
  state.names <- unlist(sapply(seq, function(seq) colnames(ecospace[[seq]]$char.space)[seq_len(ncol(ecospace[[seq]]$char.space) - 3)]))
  char.names <- sapply(seq, function(seq) ecospace[[seq]]$char)
  cs <- sapply(seq, function(seq) ncol(ecospace[[seq]]$char.space) - 3)
  for(ch in 1:nchar) {
    if(ecospace[[ch]]$type == "numeric" | ecospace[[ch]]$type == "ord.num") for(s in seq_len(cs[ch])) pre.data <- as.data.frame(cbind(pre.data, numeric(Smax)))
    if(ecospace[[ch]]$type == "factor") for(s in seq_len(cs[ch])) pre.data <- as.data.frame(cbind(pre.data, factor(Smax, ordered=FALSE, levels=ecospace[[ch]]$allowed.combos)))
    if(ecospace[[ch]]$type == "ord.fac") for(s in seq_len(cs[ch])) pre.data <- as.data.frame(cbind(pre.data, factor(Smax, ordered=TRUE, levels=ecospace[[ch]]$allowed.combos)))
  }
  pre.data <- pre.data[,2:ncol(pre.data)]
  if(length(char.names)==length(state.names)) { colnames(pre.data) <- char.names } else { colnames(pre.data) <- state.names }
  return(pre.data)
}
