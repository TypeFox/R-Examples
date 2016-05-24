#' Goodman-Kruskal's \eqn{\gamma}
#'
#' Computes Goodman-Kruskal's \eqn{\gamma}
#' @param M a matrix
#' @return
#' \item{scon}{concordance}
#' \item{sdis}{disconcordance}
#' \item{gamma}{a real number between -1 and 1. calculated as
#' \eqn{\code{gamma} = \frac{\code{scon}-\code{sdis}}{\code{scon}+\code{sdis}}}{(scon-sdis)/(scon+sdis)}}
#' @references Goodman LA, Kruskal WH (1954) Measures of association
#' for cross classifications, Journal of the American Statistical
#' Association, 49, 732-764.
#' @export

GKGamma <- function(M) {
  nrow = nrow(M)
  ncol = ncol(M)
  scon = sdis = 0
  for(i in 1:(nrow-1)) {
    for(j in 1:ncol) {
      if(j<ncol)
        scon = scon + M[i,j] * sum(M[(i+1):nrow, (j+1):ncol])
      if(j>1)
        sdis = sdis + M[i,j] * sum(M[(i+1):nrow, 1:(j-1)])
    }
  }

  list(scon=scon, sdis=sdis, gamma = (scon-sdis)/(scon+sdis))
}
