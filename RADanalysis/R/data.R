#' Abundance table of gut microbiome.
#'
#' Dethlefsen et al. 2008 (http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0060280)
#' have treated healthy individuals with the antibiotic Ciprofloxacin and monitored the states of the gut
#' microbiome before the treatment, during the treatment, and some time after the treatment.
#'
#' Abundance table of 18 gut samples from 3 individuals, prior, under and after using Ciprofloxacin (Cp)
#' antibiotic. Each row is the abundance of gut microbiome for one sample. The row names are as follows:
#'
#' \itemize{
#'   \item A1, B1 and C1 are samples from individuals A, B and C 60 days prior using Cp.
#'   \item A2a is the sample from individual A 6 days prior using Cp.
#'   \item A2b is the sample from individual A 2 days prior using Cp.
#'   \item A2c, B2 and C2 are samples from individuals A, B and C one day prior using Cp.
#'   \item A3a is the sample from individual A 3 days after the first day of Cp administration.
#'   \item A3b, B3 and C3 are the samples from individual A, B and C 5 days after the first day of Cp administration.
#'   \item A4, B4 and C4 are the samples from individual A, B and C 33 days after the first day of Cp administration.
#'   \item A5, B5 and C5 are the samples from individual A, B and C 180 days after the first day of Cp administration.
#' }
#'
#' points 1 and 2 are considered pre-Cp, points 3 are considered under-Cp and point 4 and 5 are considered post-Cp.
#'
#' @docType data
#' @keywords datasets
#' @name gut_otu_table
#' @usage data(gut_otu_table)
#' @format A matrix with 18 rows and 5670 columns.
#' @source Dethlefsen, Les, et al. "The pervasive effects of an antibiotic on the human gut microbiota, as revealed by
#'  deep 16S rRNA sequencing." PLoS biol 6.11 (2008): e280.
"gut_otu_table"


#' Result of normalization of gut_otu_table with \code{max_rank = 400} and \code{average_over = 2000}.
#'
#' Normalized rads created from \code{\link{gut_otu_table}} with \code{max_rank = 400} and \code{average_over = 2000}.
#' \code{\link{gut_otu_table}} has 18 gut samples from 3 individuals, prior, under and after using Ciprofloxacin (Cp)
#' antibiotic. Each row is the abundance of gut microbiome for one sample. The row names are as follows:
#'
#' \itemize{
#'   \item A1, B1 and C1 are samples from individuals A, B and C 60 days prior using Cp.
#'   \item A2a is the sample from individual A 6 days prior using Cp.
#'   \item A2b is the sample from individual A 2 days prior using Cp.
#'   \item A2c, B2 and C2 are samples from individuals A, B and C one day prior using Cp.
#'   \item A3a is the sample from individual A 3 days after the first day of Cp administration.
#'   \item A3b, B3 and C3 are the samples from individual A, B and C 5 days after the first day of Cp administration.
#'   \item A4, B4 and C4 are the samples from individual A, B and C 33 days after the first day of Cp administration.
#'   \item A5, B5 and C5 are the samples from individual A, B and C 180 days after the first day of Cp administration.
#' }
#'
#' points 1 and 2 are considered pre-Cp, points 3 are considered under-Cp and point 4 and 5 are considered post-Cp.
#'
#' order of rows are similar to gut_otu_table
#'
#' @docType data
#' @keywords datasets
#' @name gut_nrads
#' @usage data(gut_nrads)
#' @format Contain the result of \code{\link{RADnormalization_matrix}(input = gut_otu_table,max_rank = 400,average_over = 2000)}
"gut_nrads"
