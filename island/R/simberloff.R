#' Simberloff and Wilson original defaunation experiment data
#'
#' A list of datasets containing the presence-absence data gathered originally
#' by Simberloff and Wilson in their defaunation experiment of six mangrove
#' islands in the Florida Keys.
#'
#' The defaunation experiment of Simberloff and Wilson was aimed to test
#' experimentally the Theory of Island Biogeography.  The approach sought was
#' eliminating the fauna of several islands and following the recolonization
#' proccess.
#'
#' After some trials,  six red mangrove islets of Florida Bay were chosen for
#' the task.  These islets had to be stripped of all arthropofauna without
#' harming the vegetation and then all the colonists were identified.  The
#' result of these defaunation experiments supported the existence of species
#' equilibria and were consistent with the basic MacArthur-Wilson equilibrium
#' model.
#'
#' @note The shaded entries in the original dataset, for taxa inferred to be
#'   present from other evidence rather than direct observation, are considered
#'   as present in these datasets.
#'
#' @format A list with 6 dataframes, each corresponding to the survey of a
#'   different island. Dataframes have in columns: \describe{ \item{Taxa}{Taxa
#'   considered} \item{PRE}{Presence-absence before the defaunation process}
#'   \item{Integers (e.g. 21, 40, 58...)}{Several columns with presence-absence
#'   data for the day specified} \item{Tax. Unit 1}{Highest taxonomical unit
#'   considered} \item{Tax. Unit 2}{Second highest taxonomical unit considered}
#'   \item{Genera}{Genera of the identified taxon} \item{Island}{Island of
#'   identification of the taxon} }
#'
#' @source Simberloff, D. S., & Wilson, E. O.. (1969). Experimental Zoogeography
#'   of Islands: The Colonization of Empty Islands. \emph{Ecology},
#'   \bold{50(2)}, 278--296. \url{http://doi.org/10.2307/1934856}
#'
#' @references Wilson, E. O.. (2010). Island Biogeography in the 1960s: THEORY
#'   AND EXPERIMENT. In J. B. Losos and R. E. Ricklefs (Eds.), \emph{The Theory
#'   of Island Biogeography Revisited} (pp. 1--12). Princeton University Press.
#'   \cr \cr Simberloff, D. S., and Wilson, E. O.. (1969). Experimental
#'   Zoogeography of Islands: The Colonization of Empty Islands. \emph{Ecology},
#'   \bold{50(2)}, 278--296. \url{http://doi.org/10.2307/1934856} \cr \cr
#'   Wilson, E. O., and Simberloff, D. S.. (1969). Experimental Zoogeography of
#'   Islands: Defaunation and Monitoring Techniques. \emph{Ecology},
#'   \bold{50(2)}, 267--278. \url{http://doi.org/10.2307/1934855} \cr \cr
#'   Simberloff, D. S.. (1969). Experimental Zoogeography of Islands: A Model
#'   for Insular Colonization. \emph{Ecology}, \bold{50(2)}, 296--314.
#'   \url{http://doi.org/10.2307/1934857}
#'
#' @name simberloff
NULL
