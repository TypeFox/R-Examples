#' An overview of the function of package \code{editrules}
#'
#' 
#' The \code{editrules} package aims to provide an environment to conveniently
#' define, read and check recordwise data constraints including 
#' \itemize{
#' \item{Linear (in)equality constraints for numerical data},
#' \item{Constraints on value combinations of categorical data}
#' \item{Conditional constraints on numerical and/or mixed data}
#' }
#' In literature these constraints, or restrictions are refered to as ``edits''. 
#' \code{editrules} can perform common rule
#' set manipulations like variable elimination and value substitution, and 
#' offers error localization functionality based on the
#' (generalized) paradigm of Fellegi and Holt. Under this paradigm, one determines
#' the smallest (weighted) number of variables to adapt such that no (additional or derived) 
#' rules are violated. The paradigm is based on the assumption that errors
#' are distributed randomly over the variables and there is no detectable cause of
#' error. It also decouples the detection of corrupt variables from their
#' correction. For some types of error, such as sign flips, typing errors or
#' rounding errors, this assumption does not hold. These errors can be detected
#' and are closely related to their resolution. The reader is referred to the
#' \pkg{deducorrect} package for treating such errors. 
#'
#' @section I. Define edits:
#'
#' \code{editrules} provides several methods for creating edits from a \code{character}
#' , \code{expression}, \code{data.frame} or a text file.
#' \tabular{ll}{
#'   \code{\link{editfile}}     \tab Read  conditional numerical, numerical and categorical constraints from textfile \cr
#'   \code{\link{editset}}     \tab Create conditional numerical, numerical and categorical constraints \cr
#'   \code{\link{editmatrix}} \tab Create a linear constraint matrix for numerical data \cr
#'   \code{\link{editarray}}  \tab Create value combination constraints for categorical data \cr
#' }
#'
#' @section II. Check and find errors in data:
#'
#' \code{editrules} provides several method for checking \code{data.frame}s with edits
#' \tabular{ll}{
#'   \code{\link{violatedEdits}} \tab Find out which record violates which edit. \cr
#'   \code{\link{localizeErrors}}  \tab Localize erroneous fields using Fellegi and Holt's principle. \cr
#'   \code{\link{errorLocalizer}}  \tab Low-level error localization function using B&B algorithm \cr
#' }
#' Note that you can call \code{plot}, \code{summary} and \code{print}  on results of these functions.
#'
#' @section IV. Manipulate and check edits:
#'
#' \code{editrules} provides several methods for manipulating edits
#' \tabular{ll}{
#'   \code{\link{substValue}} \tab Substitute a value in a set of rules \cr
#'   \code{\link{eliminate}} \tab Derive implied rules by variable elimination \cr
#'   \code{\link{reduce}} \tab Remove unconstraint variables \cr
#'   \code{\link{isFeasible}} \tab Check for contradictions \cr
#'   \code{\link{duplicated}} \tab Find duplicated rules \cr
#'   \code{\link{blocks}} \tab Decompose rules into independent blocks \cr
#'   \code{\link{disjunct}} \tab Decouple conditional edits into disjunct edit sets\cr
#'   \code{\link{separate}} \tab Decompose rules in blocks and decouple conditinal edits \cr
#'   \code{\link{generateEdits}} \tab Generate all nonredundant implicit edits (\code{\link{editarray}} only) \cr
#' }
#'
#' @section V. Plot and coerce edits:
#'
#' \code{editrules} provides several methods for plotting and coercion.
#' \tabular{ll}{
#'   \code{\link{editrules.plotting}} \tab Plot edit-variable connectivity graph \cr
#'   \code{\link{as.igraph}} \tab Coerce to edit-variable connectivity \code{igraph} object \cr
#'   \code{as.character} \tab Coerce edits to \code{character} representation \cr
#'   \code{as.data.frame} \tab Store \code{character} representation in \code{data.frame} \cr
#' }
#'
#' @import lpSolveAPI
#' @importFrom igraph as.igraph
#' @name editrules-package 
#' @docType package 
{}


# not on imports: we need to DEPEND on igraph since are exporting methods of an igraph
# generic (as.igraph).


