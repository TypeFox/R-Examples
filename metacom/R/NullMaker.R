#' Null matrix creator
#'
#' Creates null matrices based on the constraints of the null model algorithm
#' ('method'). Also allows for null matrices with a species that occurs at no
#' sites, or a site without any species to be removed from the suite of
#' simulated null matrices. This function borrows heavily from the
#' commsimulator() function in the 'vegan' package, but also allows for the
#' fixed-fixed null model.
#'
#' 'method' is the null model algorithm used to create the null matrices. The
#' choice of a null algorithm is nontrivial. Leibold & Mikkelson advocated the
#' use of equiprobable rows and columns (provided that rows and columns had at
#' least one entry). This method is called 'r00'. Methods maintaining row
#' (site) frequencies include 'r0','r1' & 'r2', whereas species (column)
#' occurrences are preserved with fixed column methods such as 'c0'. The
#' default method argument is 'r1', which maintains the species richness of a
#' site (row totals) and fills species ranges (columns) based on their marginal
#' probabilities. Arguably the most conservative null algorithm is the fixed
#' row - fixed column total null, which can be attained using many of swap
#' algorithms described in the vegan package (sequential methods like 'tswap',
#' 'swap', and non-sequential 'quasiswap' and 'backtracking'). Other
#' randomization methods are also available. See the help file for 'commsim',
#' or Wright et al. 1998 for more information.
#'
#' @param comm community data in the form of a presence absence matrix
#' @param sims number of simulated null matrices to use in analysis
#' @param method null model randomization method. See details below.
#' @param ordinate logical. Would you like to ordinate the null matrices?
#' Default is TRUE.
#' @param scores Axis scores to ordinate matrix. 1: primary axis scores
#' (default) 2: secondary axis scores. See Details.
#' @param allowEmpty logical argument indicating whether to allow null
#' matrices to have empty rows or columns
#' @param verbose Logical. Prints a graphical progress bar that tracks the
#' creation of null matrices. Useful for conservative null models on large
#' and/or sparse data.
#' @return rmats -- A list of length(sim) containing the null matrices
#' @author Tad Dallas and John Lefcheck
#' @export
#' @seealso commsimulator(), nullmodel(), permatfull(), commsim()
#' @references J. Oksanen, F.G. Blanchet, R. Kindt, P. Legendre, P.R. Minchin,
#' R.B. O'Hara, G.L. Simpson, P. Solymos, M.H.H. Stevens and H. Wagner (2012).
#' vegan: Community Ecology Package. R package version 2.0-4.
#' http://CRAN.R-project.org/package=vegan
#' @keywords ordination
#' @examples
#'
#' #define an interaction matrix
#' data(TestMatrices)
#' intmat <- TestMatrices[[7]]
#'
#' #creation of the null matrices
#' nulls <- NullMaker(intmat, sims=100, method='r1')
#'

NullMaker = function (comm, sims = 1000, method = "r1", ordinate = TRUE, scores = 1, allowEmpty = FALSE, verbose = FALSE) {

  if(verbose == TRUE) pb = txtProgressBar(min = 0, max = sims, style = 3)

  #generate null matrices
  nm <- nullmodel(comm, method = method)
  sm <- simulate(nm, nsim = sims)
  sm.list <- lapply(seq(dim(sm)[3]), function(i) sm[, , i])

  if(verbose == TRUE & allowEmpty == TRUE) setTxtProgressBar(pb, length(sm.list))

  if(allowEmpty == FALSE) {
    flag = FALSE
    #remove matrices with empty rows or cols
    sm.list <- if(allowEmpty == FALSE) lapply(sm.list, function(i) if(any(colSums(i) == 0) | any(rowSums(i) == 0)) NULL else i)
    sm.list[sapply(sm.list,is.null)] = NULL

    if(verbose == TRUE)  setTxtProgressBar(pb, length(sm.list))
    while(flag == FALSE) {
      if(length(sm.list) == sims) flag = TRUE else {
        #generate extra matrices
        spares <- simulate(nm, nsim = sims/10)
        spares.list <- lapply(seq(dim(spares)[3]), function(i) spares[, , i])
        spares.list <- lapply(spares.list, function(i) if(any(colSums(i) == 0) | any(rowSums(i) == 0)) NULL else i)
        spares.list[sapply(spares.list,is.null)] <- NULL

        #replace in original sm object
        sm.list <- append(sm.list, spares.list)

        if(verbose == TRUE & length(sm.list) <= sims)  setTxtProgressBar(pb, length(sm.list))

        if(length(sm.list) >= sims) {
          sm.list <- sm.list[1:sims]
          if(verbose == TRUE)  setTxtProgressBar(pb, length(sm.list))
          flag = TRUE }
      }
    }
  }

  #Run ordination
  if(ordinate == TRUE) sm.list <- lapply(sm.list, OrderMatrix, scores = scores)
  return(sm.list)
}
