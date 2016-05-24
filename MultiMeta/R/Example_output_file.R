#' Example output file
#'
#'@description This file is an example of the output obtained with \code{multi_meta} function.
#'The input files were Example_file_1.txt and Example_file_2.txt, the command line was the one in the example (see \code{multi_meta} help)
#'
#'The columns are:
#'\itemize{
#'\item{\strong{chr} Chromosome}
#'\item{\strong{SNP} SNP name} 
#'\item{\strong{Position} Position}
#'\item{\strong{allele1} Effect allele} 
#'\item{\strong{allele0} Non-effect allele} 
#'\item{\strong{tot_af} Effect-allele frequency computed on all the cohorts analysed}
#'\item{\strong{n_pops} Number of cohorts (populations) analysed}
#'\item{\strong{pops} String of 1 or 0 values indicating which populations were analysed: e.g. if the first number is 1, the first
#'cohort was included in the analysis etc.} 
#'\item{\strong{beta_1}, \strong{beta_2}, ..., \strong{beta_p} Effect sizes for each of the \emph{p} traits; in this example \emph{p=6} }
#'\item{\strong{Vbeta_1_1}, \strong{Vbeta_1_2}, ..., \strong{Vbeta_1_p}, \strong{Vbeta_2_2}, ...,
#'\strong{Vbeta_2_p}, ..., \strong{Vbeta_p_p} Variance-covariance matrix entries (diagonal and upper triangle values only, since this matrix is symmetric)}
#'\item{\strong{p_value} P-value}
#'}
#'
#' @name Example_output_file
#' @docType data
#' @author Dragana Vuckovic \email{dragana.vuckovic@@burlo.trieste.it}
#' @keywords data
NULL     