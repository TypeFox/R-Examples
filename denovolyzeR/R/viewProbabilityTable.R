#' Displays underlying \emph{de novo} probability tables
#'
#' Tabulates probability of \emph{de novo} variant for each protein-coding variant class, for each gene.  Values are probability of a \emph{de novo} variant per chromosome per generation.  i.e. expected number of de novos for a given gene/class = \eqn{p * 2 * nsamples}.
#' @param format option to display table in wide format (default; one line per gene), or long format
#' @export

# --------------------

viewProbabilityTable <- function(format="wide"){
  if(format=="long"){
    return(pDNM)
  } else if (format=="wide") {
    return(
      reshape::cast(pDNM, hgncID + hgncSymbol + enstID + ensgID + geneName ~ class)
      )
  }
}
