#' Generate pedigrees from multi-parent designs
#' 
#' Generates pedigrees for 4-way or 8-way MAGIC designs, with specifications given by the user. These pedigrees can be used to simulate datasets with other functions in the mpMap package.
#' @export sim.mpped
#' @aliases gen4ped gen8ped
#' @param nfounders Number of founders in the pedgiree
#' @param nfunnels Number of different permutations of founder crosses (funnels) to generate.  
#' @param nperfam Number of plants per funnel (before selfing begins)
#' @param nssdgen Number of generations of single-seed descent (selfing)
#' @param nseeds Number of seeds propagated from each plant through SSD
#' @param iripgen Number of generations of advanced intercrossing between mixing and selfing
#' @param seed Random seed for selecting funnels to generate
#' @param \dots Additional arguments
#' @details The 4-way pedigree is relatively simple compared to the 8-way pedigree, as there are many fewer possible combinations of lines to intercross. For a 4-way cross there is a maximum of 3 possible funnels (ABCD, ACBD and ADBC), while there are 315 for an 8-way cross. This is computed under the assumption that cross AB is equivalent to cross BA. 
#'
#' For a value of 1, the funnel will be of the form ABCD or ABCDEFGH. For values less than the maximum possible number of funnels the crosses will be selected randomly, with each possibly funnel having equal probability of being generated. 
#' @return Matrix with four columns: ID, mother, father, and observed status
#' @seealso \code{\link[mpMap]{sim.mpcross}}
#' @examples
#' ped4way <- sim.mpped(4, 1, 500, 6, 1)
#' ped8way <- sim.mpped(8, 1, 500, 6, 1)

sim.mpped <-
function(nfounders, nfunnels=1, nperfam=50, nssdgen=6, nseeds=1, iripgen=0, seed=1, ...)
{
  if (missing(nfounders))
	stop("Must input the number of founders in the pedigree")

  if (nfounders==4)
	ped <- gen4ped(nfunnels, nperfam, nssdgen, nseeds, iripgen,...)

  if (nfounders==8)
	ped <- gen8ped(nfunnels, nperfam, nssdgen, nseeds, iripgen, seed, ...)

  attr(ped, "call") <- match.call()

  return(ped)
}

