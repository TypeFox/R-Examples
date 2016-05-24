#' Count how many generations of advanced intercross are in a pedigree
#'
#' Counts the number of generations of breeding preceding selfing and 
#' subtracts off the number necessary to minimally mix the founders' genomes
#' @export
#' @param pedigree Pedigree for a multi-parent cross. Can be generated using \code{\link[mpMap]{sim.mpped}}
#' @return Integer - number of generations of advanced intercrossing after mixing stage but before selfing.
#' @seealso \code{\link[mpMap]{sim.mpped}}
#' @examples
#' sim.map <- list(Chr1=seq(0,100,10))
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' nai(sim.ped)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1, 5)
#' nai(sim.ped)

nai <- function(pedigree)
{
  # returns the number of intercrossing generations before selfing
  nped <- nrow(pedigree)

  id <- nped
  while (pedigree[id, 2]==pedigree[id, 3])
	id <- pedigree[id,2]
 
  ngen <- 0
  while (pedigree[id,2]>0)
  {
	id <- pedigree[id,2]
	ngen <- ngen+1
  }

  nfounders <- sum(pedigree[,2]==0 & pedigree[,3]==0)
  ngen <- ngen-log(nfounders)/log(2)

  return(ngen)
}
