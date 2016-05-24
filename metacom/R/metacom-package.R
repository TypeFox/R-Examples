#' @title Pattern-based analysis of metacommunity structure
#'
#' @description   'metacom' provides functions for the analysis of the elements of
#' metacommunity structure (coherence, boundary clumping, & turnover),
#' following the pattern-based metacommunity framework of Leibold & Mikkelson
#' 2002 and Presley et al. 2010. This package is designed to allow the user to
#' distinguish between several idealized patterns of metacommunity structure
#' (Presley et al. 2010) utilizing any number of null model algorithms for the
#' randomization procedure. However, these metrics can also be used in
#' isolation, and without ordination via reciprocal averaging, and instead,
#' ordering along some biological gradient.
#'
#'
#' A metacommunity is a set of sites (e.g. plants in plant-pollinator networks)
#' associated through interactions (e.g. insect species (columns) interact with
#' plant species (rows) in plant-pollinator networks). The pattern-based
#' metacommunity concept, proposed by Leibold & Mikkelson 2002 and expounded on
#' by Presley et al. 2010, allows for the evaluation of metacommunity structure
#' by using randomization techniques to discern between 10 patterns of
#' metacommunity structure. This is performed by ordinating site-by-species
#' interaction matrices and calculating three metrics; coherence, boundary
#' clumping & turnover.
#'
#' The metacom package calculates these three metrics; coherence is calculated
#' using the function Coherence(), boundary clumping with BoundaryClump(), and
#' turnover (from either species or range perspective) using the Turnover()
#' function. These functions are consolidated in the Metacommunity() function,
#' which can be used to calculate all three metrics. In order to interpret the
#' output of these functions, it will be helpful to read Leibold & Mikkelson
#' 2002 and Presley et al. 2010, but to also read Ulrich and Gotelli 2013, as
#' this paper outlines the difficulty seemingly inherent with investigating
#' community structure. Also, these functions do not have to be used strictly
#' in the Leibold and Mikkelson 2002 framework.
#'
#' I caution the user to be aware that the creation of null matrices can be
#' performed to allow (or not allow) sites to be empty, or species to not exist
#' at any site (i.e. column sums and/or row sums are allowed to be zero). This
#' is controlled by the logical argument 'allow.empty' in the Metacommunity(),
#' NullMaker(), Coherence(), and Turnover() functions. Restricting nulls to not
#' allow empty rows or columns may be biologically realistic, but it also
#' reduces the number of unique null matrices that can be built, which will
#' impact computation time, making it infeasible or impossible in some
#' situations. These situations occur when you have a very sparse interaction
#' matrix, and is also influenced by null model algorithm ('method') that you
#' choose.
#'
#' The 'metacom' package is partially adapted from previous Matlab code written
#' by Christopher Higgins (available at
#' http://www.tarleton.edu/Faculty/higgins/EMS.htm) and relies on many
#' functions in the 'vegan' package (Oksanen et al. 2012)
#'
#' @name metacom-package
#' @aliases metacom-package metacom
#' @docType package
#' @author Tad Dallas
#' @import vegan
#' @importFrom graphics axis box image mtext par
#' @importFrom stats pchisq pnorm sd simulate
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @references
#'
#' Dallas,T. 2014. metacom: an R package for the analysis of metacommunity
#' structure. Ecography. DOI:10.1111/j.1600-0587.2013.00695.x
#'
#' Leibold, M. & Mikkelson, G. (2002) Coherence, species turnover, and boundary
#' clumping: elements of metacommunity structure. Oikos, 97, 237-250.
#'
#' Leibold, M., Holyoak, M., Mouquet, N., Amarasekare, P., Chase, J., Hoopes,
#' M., Holt, R., Shurin, J., Law, R., Tilman, D. et al. (2004) The
#' metacommunity concept: a framework for multi-scale community ecology.
#' Ecology letters, 7, 601-613.
#'
#' Oksanen, J., F.G. Blanchet, R. Kindt, P. Legendre, P.R. Minchin, R.B.
#' O'Hara, G.L. Simpson, P. Solymos, M.H.H. Stevens and H. Wagner (2012).
#' vegan: Community Ecology Package. R package version 2.0-4.
#' http://CRAN.R-project.org/package=vegan
#'
#' Presley, S., Higgins, C. & Willig, M. (2010) A comprehensive framework for
#' the evaluation of metacommunity structure. Oikos, 119, 908-917.
#'
#' Ulrich, W. and Gotelli, N. J. (2013) Pattern detection in null model
#' analysis. Oikos, 122: 2-18. doi: 10.1111/j.1600-0706.2012.20325.x
#'
#' Willig, M., Presley, S., Bloch, C., Castro-Arellano, I., Cisneros, L.,
#' Higgins, C. & Klingbeil, B. (2011) Tropical metacommunities along
#' elevational gradients: effects of forest type and other environmental
#' factors. Oikos, 120, 1497-1508.
NULL





#' Test matrices used to evaluate metacommunity functions
#'
#' A list of 7 test matrices from two of the foundational papers on the
#' Elements of Metacommunity Structure analysis (Leibold & Mikkelson 2002 and
#' Presley et al. 2010)
#'
#'
#' @name TestMatrices
#' @docType data
#' @format A list containing interaction matrices from Leibold & Mikkelson 2002
#' and Presley et al. 2010:
#'
#' 1) dim=(20 x 20) Randomly generated matrix (`rbinom(400,1,0.4)`)
#'
#' 2) dim=(10 x 10) Leibold & Mikkelson 2002 Figure 1b
#'
#' 3) dim=(10 x 10) Leibold & Mikkelson 2002 Figure 2a
#'
#' 4) dim=(10 x 10) Leibold & Mikkelson 2002 Figure 2b
#'
#' 5) dim=(15 x 10) Leibold & Mikkelson 2002 Figure 3c
#'
#' 6) dim=(20 x 20) Presley et al. Figure 3c
#'
#' 7) dim(20 x 20) Presley et al. Figure 4a
#'
#' @source Leibold, M. A., & Mikkelson, G. M. (2002). Coherence, species
#' turnover, and boundary clumping: elements of metacommunity structure. Oikos,
#' 97(2), 237-250.
#'
#' Presley, S. J., C. L. Higgins, and M. R. Willig. 2010. A comprehensive
#' framework for the evaluation of metacommunity structure. Oikos 119:908-917
#' @keywords datasets
#' @examples
#'
#' #load list containing interaction matrices
#' data(TestMatrices)
#'
#' length(TestMatrices)
#' names(TestMatrices)
#'
#' #image plot of interaction matrix, using the Imagine() function
#' test <- TestMatrices[[6]]
#' Imagine(test, fill=FALSE, order=TRUE)
#'
#'
NULL
