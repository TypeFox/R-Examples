#' palinear2LiLog
#'
#' Apply LiLog to \code{palinear} output files
#'
#' The function tranforms betas from the \code{palinear} output files
#' into their logistic correspondents.
#'
#' @param prevalence: prevalence of the studied diseases
#' @param input.file: palinear output file that needs transformation
#' @param output.file: name of the transformed file, if NULL
#'  \code{input.file} will be overwritten
#'
#' @return Transformed file
#'
#' @author Nicola Pirastu, Lennart Karssen
#'
#' @seealso
#' For futher details on the methods and for the basic function see
#' \code{\link{LiLog}}. For transforming the output from
#' \code{formetascore} use \code{\link{formetascore2LiLog}}.
#'
#'
#' @examples
#'
#' \dontrun{
#'  palinear2LiLog(prevalence=0.2, input.file="chr18.add.out.txt",
#'                 output.file="chr18.logistic.add.out.txt")
#' }
#'
#' @references
#'
#' \emph{Genome-wide feasible and meta-analysis oriented methods for
#' analysis of binary traits using mixed models},
#' Nicola Pirastu, Lennart C. Karssen, Pio D'adamo, Paolo Gasparini,
#' Yurii Aulchenko (Submitted).


palinear2LiLog <- function(prevalence, input.file, output.file=NULL){
  if(is.null(output.file)){
    output.file <- input.file
  }

  res  <- read.table(input.file, header=TRUE, stringsAsFactors=FALSE)
  Lres <- LiLog(prev=prevalence, beta=res$beta_SNP_add,
                SEbeta=res$sebeta_SNP_add, q=res$Mean_predictor_allele/2.)

  res$beta_SNP_add   <- Lres$beta
  res$sebeta_SNP_add <- Lres$sebeta

  write.table(res, file=output.file, row.names=FALSE, quote=FALSE)
}
