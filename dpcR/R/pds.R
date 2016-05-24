#' Plasmid dilution series results
#' 
#' These are the results data from the \code{pds_raw} data as calculated by the
#' BioRad QX100 Droplet Digital PCR System.
#' 
#' Setup: Duplex assay with constant amount of genomic DNA and six 10-fold
#' dilutions of plasmid DNA with 4 replicates, ranging theoretically from ~
#' 10^4 to 10^-1 copies/ micro L plus 4 replicates without plasmid DNA.
#' Included are No-gDNA-control and No-template-control, 2 replicates each.
#' 
#' Annotation: FX.Y (X = dilution number, Y = replicate number). Hardware:
#' Bio-Rad QX100 Droplet digital PCR system Details: Genomic DNA isolated from
#' Pseudomonas putida KT2440. Plasmid is pCOM10-StyA::EGFP StyB [Jahn et al.,
#' 2013, Curr Opin Biotechnol, Vol. 24 (1): 79-87]. Template DNA was heat
#' treated at 95 degree Celsius for 5 min prior to PCR. Channel 1, primers for
#' genomic DNA marker ileS, Taqman probes (FAM labelled). Channel 2, primers
#' for plasmid DNA marker styA, Taqman probes (HEX labelled).
#' 
#' 
#' @name pds
#' @docType data
#' @format A data frame with 64 observations on the following 44 variables.
#' \describe{ 
#' \item{Well}{a factor with levels \code{A01} to \code{H04}} 
#' \item{ExptType}{a factor with levels \code{Absolute Quantification}} 
#' \item{Experiment}{a factor with levels \code{ABS}}
#' \item{Sample}{a factor with levels \code{B} \code{B + P 10^2}
#' \code{gDNA} \code{gDNA + P 10^0} \code{gDNA + P 10^1} \code{gDNA + P 10^-1}
#' \code{gDNA + P 10^2} \code{gDNA + P 10^3} \code{gDNA + P 10^4}}
#' \item{TypeAssay}{a factor with levels \code{Ch1NTC}
#' \code{Ch1Unknown} \code{Ch2NTC} \code{Ch2Unknown}} 
#' \item{Assay}{a factor with levels \code{ileS} \code{styA}} 
#' \item{Status}{a factor with levels \code{Manual}} 
#' \item{Concentration}{a numeric vector}
#' \item{TotalConfMax}{a logical vector} 
#' \item{TotalConfMin}{a logical vector} 
#' \item{PoissonConfMax}{a numeric vector}
#' \item{PoissonConfMin}{a numeric vector} 
#' \item{Positives}{a numeric vector} 
#' \item{Negatives}{a numeric vector}
#' \item{Ch1.Ch2.}{a numeric vector} 
#' \item{Ch1.Ch2..1}{a numeric vector} 
#' \item{Ch1.Ch2..2}{a numeric vector}
#' \item{Ch1.Ch2..3}{a numeric vector} 
#' \item{Linkage}{a numeric vector} 
#' \item{AcceptedDroplets}{a numeric vector}
#' \item{CNV}{a logical vector} 
#' \item{TotalCNVMax}{a logical vector} 
#' \item{TotalCNVMin}{a logical vector}
#' \item{PoissonCNVMax}{a logical vector}
#' \item{PoissonCNVMin}{a logical vector}
#' \item{ReferenceCopies}{a logical vector}
#' \item{UnknownCopies}{a logical vector} 
#' \item{Ratio}{a numeric vector} 
#' \item{TotalRatioMax}{a logical vector}
#' \item{TotalRatioMin}{a logical vector}
#' \item{PoissonRatioMax}{a numeric vector}
#' \item{PoissonRatioMin}{a numeric vector}
#' \item{FractionalAbundance}{a numeric vector}
#' \item{TotalFractionalAbundanceMax}{a logical vector}
#' \item{TotalFractionalAbundanceMin}{a logical vector}
#' \item{PoissonFractionalAbundanceMax}{a numeric vector}
#' \item{PoissonFractionalAbundanceMin}{a numeric vector}
#' \item{ReferenceAssayNumber}{a numeric vector}
#' \item{TargetAssayNumber}{a numeric vector}
#' \item{MeanAmplitudeofPositives}{a numeric vector}
#' \item{MeanAmplitudeofNegatives}{a numeric vector}
#' \item{MeanAmplitudeTotal}{a numeric vector}
#' \item{ExperimentComments}{a logical vector}
#' \item{MergedWells}{a logical vector} }
#' @author Michael Jahn, Stefan Roediger, Michal Burdukiewcz
#' @references Jahn et al., 2013, \emph{Curr Opin Biotechnol}, Vol. 24 (1):
#' 79-87
#' @source Michael Jahn Flow cytometry group / Environmental microbiology
#' Helmholtz Centre for Environmental Research - UFZ Permoserstrasse 15 / 04318
#' Leipzig / Germany phone +49 341 235 1318 michael.jahn [at] ufz.de /
#' www.ufz.de
#' @keywords datasets
#' @examples
#' summary(extract_dpcr(read_QX100(pds), 1L:10))
#' 
NULL