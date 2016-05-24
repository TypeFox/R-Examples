#' Coverage dataset for FLO8 gene.
#'
#' A dataset containing the RNA-seq coverage values
#' in all samples for FLO8 gene.
#'
#' @format An object of type geneCoverage.
#'
"yer109c"

#' Coverage dataset for FLO8, CDC14 and SNC1 gene.
#'
#' A dataset containing the RNA-seq coverage values
#' in all samples for FLO8, CDC14 and SNC1 gene.
#'
#' @format An object of type coverageDataSet.
#'
"rna_seq_signals"

#' Saccharomyce cerevisiae genes annotations.
#'
#' A dataset containing a subset of four genes from saccer3 gene annotations.
#'
#' @format A data.table with 4 rows and 5 variables.
#'
#' \describe{
#'   \item{chromosome}{Gene's chromosome name}
#'   \item{start}{Gene's transcript start site}
#'   \item{end}{Gene's transcript termination site}
#'   \item{gene_name}{Gene's SGD name}
#'   \item{strand}{Gene's strand}
#' }
#'
"annotations"

#' Saccharomyce cerevisiae introns annotations.
#'
#' A dataset containing saccer3 introns annotations.
#'
#' @format A data.table with 337 rows and 5 variables.
#'
#' \describe{
#'   \item{chromosome}{Intron' chromosome name}
#'   \item{start}{Intron's transcript start site}
#'   \item{end}{Intron's transcript termination site}
#'   \item{gene_name}{Gene's SGD name}
#'   \item{strand}{Intron's strand}
#' }
#'
"introns"

