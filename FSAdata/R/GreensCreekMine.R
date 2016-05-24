#' @title Catches in removal events of Coho Salmon and Dolly Varden Char at various locations near the Greens Creek (AK) Mine site.
#' 
#' @description Catches in removal events of Coho Salmon (\emph{Oncorhynchus kisutch}) and Dolly Varden Char (\emph{Salvelinus malma}) at various locations near the Greens Creek (AK) Mine site.
#' 
#' @details Reaches were isolated by natural features, such as shallow riffles.  The sample reaches were saturated with 6.35 mm (0.25 in) minnow traps baited with whirl packs containing disinfected salmon eggs.  The traps were deployed for 1.5 h and then retrieved where each fish was transferred into a plastic bucket, and the trap was re-baited ans re-set for another 1.5 h soak.  In between trapping events, fish were processed -- measured and recorded FL to the nearest 1 mm, weight to the nearest 0.1 g, and species identified. Captured fish were retained during the sample period and returned alive after all three passes were complete.
#' 
#' @name GreensCreekMine
#' 
#' @docType data
#' 
#' @format A data frame of 66 observations on the following 8 variables:
#'  \describe{
#'    \item{location}{Sampling location.}
#'    \item{year}{Sampling year.}
#'    \item{species}{Species (\code{Coho.Salmon} or \code{Dolly.Varden}).}
#'    \item{set1}{Catch on the first removal pass.}
#'    \item{set2}{Catch on the second removal pass.}
#'    \item{set3}{Catch on the third removal pass.}
#'    \item{min.FL}{Minimum observed fork length.}
#'    \item{max.FL}{Maximum observed fork length.}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Population size
#'    \item Abundance
#'    \item Removal
#'  }
#'  
#' @concept Abundance 'Population Size' Removal
#' 
#' @source From Appendix C1 of Kanouse, K.M. and B.P. Brewster.  2012.  Aquatic Biomonitoring at Greens Creek Mine, 2012.  Alaska Department of Fish and Game Technical Report No. 12-11.  [Was (is?) from http://dnr.alaska.gov/mlw/mining/largemine/greenscreek/pdf/gc2012bio.pdf.]
#' 
#' @keywords datasets
#' 
#' @examples
#' data(GreensCreekMine)
#' str(GreensCreekMine)
#' head(GreensCreekMine)
#' 
#' ## extract data for one location, year, and species (e.g., 3rd row)
#' GreensCreekMine[3,]
#' 
NULL