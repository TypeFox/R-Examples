#' qiimer: Read QIIME output files and create plots
#' 
#' @name qiimer
#' @docType package
#' @import grDevices
#' @import pheatmap
#' @import stats
#' @import utils
NULL


#' Sample dataset from murine gut microbiome
#' 
#' The \code{relmbeta} dataset is taken from a mouse study where wild-type 
#' and RELMbeta knockout mice were fed either a normal or high-fat diet.  
#' The diet was observed to have a pronounced effect on the gut microbiome 
#' composition.  The genotype also had an effect, but less so.
#' 
#' The bacterial 16S rDNA gene was sequenced using 454 FLX technology, 
#' producing about 26k reads.  The reads were processed via the standard QIIME
#' workflow, \code{pick_de_novo_otus.py}, which clustered the reads into 776 
#' operational taxonomic units (OTUs).  QIIME version 1.7.0 was used for the 
#' analysis.
#' 
#' The OTU table was filtered to remove OTUs appearing in only one sample.
#' Following this, a single rarefaction was performed at a level of 500 reads
#' per sample.  The unweighted UniFrac distance was then computed for each pair 
#' of samples.
#' 
#' The \code{relmbeta} data frame lists the sample IDs, genotypes, and dietary 
#' assignments for the mice.
#' 
#' @name relmbeta
#' @docType data
#' @keywords datasets
#' @format A data frame with 20 rows and 3 columns.
#' @source Hildebrandt et al. High-fat diet determines the composition of the 
#'   murine gut microbiome independently of obesity. \emph{Gastroenterology} 
#'   \strong{137}, 1716 (2009).
NULL


#' OTU counts from murine gut microbiome
#' 
#' The matrix \code{relmbeta_counts} contains the number of reads 
#' observed in each OTU, after rarefaction to 500 reads per sample.  
#' OTUs are listed in the rows and samples are listed in the columns.
#'
#' The \code{\link{relmbeta}} documentation provides an overview of the study.
#' 
#' @name relmbeta_counts
#' @docType data
#' @keywords datasets
#' @format An integer matrix with 337 rows and 20 columns.
#' @seealso \code{\link{relmbeta}}
#' @source Hildebrandt et al. High-fat diet determines the composition of the 
#'   murine gut microbiome independently of obesity. \emph{Gastroenterology} 
#'   \strong{137}, 1716 (2009).
NULL


#' Unweighted UniFrac distances from murine gut microbiome
#' 
#' \code{relmbeta_dist} is an object of class \code{"dist"}, containing
#' the unweighted UniFrac distances between samples.
#'
#' The \code{\link{relmbeta}} documentation provides an overview of the study.
#' 
#' @name relmbeta_dist
#' @docType data
#' @keywords datasets
#' @format An object of class \code{"dist"}.
#' @seealso \code{\link{relmbeta}}
#' @source Hildebrandt et al. High-fat diet determines the composition of the 
#'   murine gut microbiome independently of obesity. \emph{Gastroenterology} 
#'   \strong{137}, 1716 (2009).
NULL


#' Shannon diversity measurements from murine gut microbiome
#' 
#' \code{relmbeta_alpha} is a data frame containing Shannon diversity 
#' (base 2) at a level of 10-500 sequences per sample.
#'
#' The \code{\link{relmbeta}} documentation provides an overview of the study.
#' 
#' @name relmbeta_alpha
#' @docType data
#' @keywords datasets
#' @format A data frame with 2200 rows and 4 columns.
#' @seealso \code{\link{relmbeta}}
#' @source Hildebrandt et al. High-fat diet determines the composition of the 
#'   murine gut microbiome independently of obesity. \emph{Gastroenterology} 
#'   \strong{137}, 1716 (2009).
NULL


#' Taxonomic assignments from murine gut microbiome
#' 
#' The character vector \code{relmbeta_assignments} contains taxonomic
#' assignments for each OTU in the study.
#'
#' The \code{\link{relmbeta}} documentation provides an overview of the study.
#' 
#' @name relmbeta_assignments
#' @docType data
#' @keywords datasets
#' @format A character vector with 337 elements.
#' @seealso \code{\link{relmbeta}}
#' @source Hildebrandt et al. High-fat diet determines the composition of the 
#'   murine gut microbiome independently of obesity. \emph{Gastroenterology} 
#'   \strong{137}, 1716 (2009).
NULL


#' BIOM format object from murine gut microbiome
#' 
#' The \code{relmbeta_biom} object is a representation of the BIOM file
#' produced by QIIME.  It was produced by loading the BIOM file with the 
#' function \code{RJSONIO::fromJSON}.
#'
#' The \code{\link{relmbeta}} documentation provides an overview of the study.
#' 
#' @name relmbeta_biom
#' @docType data
#' @keywords datasets
#' @format A \code{"biom"} object with 227 rows and 20 columns.
#' @seealso \code{\link{relmbeta}}
#' @source Hildebrandt et al. High-fat diet determines the composition of the 
#'   murine gut microbiome independently of obesity. \emph{Gastroenterology} 
#'   \strong{137}, 1716 (2009).
NULL
