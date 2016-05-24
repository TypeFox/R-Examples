
#' Create a geneAnnotation object.
#'
#' Create a geneAnnotation object containing the informations on
#' chromosomes, start and end positions, strand, score and name of genes.
#'
#' @param biomart	 BioMart database name you want to connect to.
#' @param dataset	Dataset you want to use.
#' @param host Host to connect to.
#'
#'
#' @details Please see \code{\link[biomaRt]{useMart}} for further information.
#'
#' @export
#' @return An object of class geneAnnotation with the following components:
#'        \item{chromosome}{Gene chromosome.}
#'        \item{start}{Gene transcription start site.}
#'        \item{end}{Gene transcription termination site.}
#'        \item{name}{Gene name.}
#'        \item{score}{Gene score.}
#'        \item{strand}{Gene strand.}
#'
#' @examples
#' annotationsSet()
annotationsSet <- function(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "scerevisiae_gene_ensembl",
                           host = "www.ensembl.org"){
  mart <- biomaRt::useMart(biomart = biomart, dataset = dataset,
                           host = host)
  annotations <- biomaRt::getBM(
    attributes = c("chromosome_name", "start_position",
                   "end_position","ensembl_gene_id", "strand"),
    mart = mart
  )
  annotations <- data.table::data.table(
    annotations[1:4], score = ".", annotations[5])
  annotations[, strand := ifelse(strand == "-1", "-", "+")]
  as.annotationsSet(annotations)
}

#' Convert a data.table into a geneAnnotation object.
#'
#' Create a geneAnnotation object containing the informations on
#' chromosomes, start and end positions, strand, score and name of genes.
#'
#' @param data A \code{\link[data.table]{data.table}} with 6 columns corresponding to the genes'
#' chromosome, start posittion, end position, gene name, gene score and gene strand.
#'
#' @export
#' @return An object of class geneAnnotation with the following components:
#'        \item{chromosome}{Gene chromosome.}
#'        \item{start}{Gene transcription start site.}
#'        \item{end}{Gene transcription termination site.}
#'        \item{name}{Gene name.}
#'        \item{score}{Gene score.}
#'        \item{strand}{Gene strand.}
#'
#' @examples
#' data(annotations)
#' as.annotationsSet(annotations)
as.annotationsSet <- function(data) {
  if (!data.table::is.data.table(data)) {
    stop("You must provide the annotations in a data.table")
  }
  if (length(data) != 6) {
    stop(
      "The data must be in a BED format. It must contains six variables:
      chromosome, start_position, end_position, gene_name, score, strand."
    )
  }

  data.table::setnames(data, c("chromosome", "start", "end", "name", "score", "strand"))

  if (!(unique(data[,strand]) == c("-", "+") ||
        unique(data[,strand]) == c("+", "-"))) {
    stop ("The last column doesn't define the strand. Must be either \"+\" or \"-\".")
  }

  data.table::setkey(data, name)
  class(data) <- append(class(data), "annotationsSet")
  data
}

#select a specified gene name from a geneAnnotations object.
.get_annotation <- function(annotationsSet, name){
  if(!is.element( "annotationsSet", class(annotationsSet)))
    annotationsSet <- as.annotationsSet(annotationsSet)
  annotationsSet[name]
}

.get_genes<- function(annotationsSet){
  if(!is.element( "annotationsSet", class(annotationsSet)))
    annotationsSet <- as.annotationsSet(annotationsSet)
  annotationsSet[,name]
}
