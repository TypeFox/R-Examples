# !/usr/bin/env R

# Version 1.0
# Author Nicole Uwimana
# Date: 02/02/2015
# Last modified 23/10/2015


# write data into file
#' @importFrom utils write.table
.writeFile <- function(data, file, append = FALSE, col.names = FALSE) {
  write.table(
    data, file, quote = F, append = append,
    sep = "\t", col.names = col.names, row.names = F
  )
}

#   string manipulation --------------------------------------------------------------------

# substr_left a character vector from the begining.
#
# substr_left a character vector from the start to the nth character.
#
# @param string a character vector.
# @param n an integer indicating the number of characters to substr_left.
#          if negative value, the string is returned without the last \code{n} characters.
#
# @return The first \code{n} characters of the \code{string} or
# the \code{string} without the last \code{n} characters.
#
# @examples
# yCrypticRNAs:::substr_left("output.txt", -4)
# yCrypticRNAs:::substr_left("output.txt", 6)
substr_left <- function(string, n){
  n <- as.numeric(n)
  if (n < 0) {
    base::substr(string, start = 1, stop = base::nchar(string) + n)
  } else {
    base::substr(string, start = 1, stop = n)
  }
}

#   vector, data.frame and data.table manipulations -----------------------------------------

# Scale data.
#
#
# @param data a data.frame, data.table or matrix of values to scale or a input file name.
# @param column an integer indicating the column of the data to scale.
# @param sf a number indicating the scaling factor.
#
# @return scaled \code{data}
#
# @import data.table
#
# @examples
# require(data.table)
# data("annotations")
# annotations[, score := sample(10000, nrow(annotations))]
# annotations
# yCrypticRNAs:::scale_data(data = annotations, column = 5, sf = 0.02)
# annotations
scale_data <- function(data, column = 8, sf){

  if(is.character(data)){
    file <- data
    data <- tryCatch(
      data.table::fread(file),
      error = function(e) stop("From scale_data function: ", e)
    )
  }else{
    data <- data.table::as.data.table(data)
    file <- NULL
  }

  data[, (column) := data[,column, with = F] * sf]

  if(!is.null(file)){
    .writeFile(data, file)
    return(paste("see results in file: ", file))
  }
  data
}

#   bedtools  coverage ------------------------------------------------------------------------

#' Report coverage.
#'
#' Report the coverage using \href{http://bedtools.readthedocs.org/en/latest/content/tools/coverage.html}{bedtools coverage} algorithm.
#'
#' @param reads a \code{\link[data.table]{data.table}} of values
#'        or a character vector indicating the input file name.
#'        Note: the data must be in a BED-like format.
#' @param annotation a object of type \code{\link{annotationsSet}}
#'        or a character vector indicating the input BED file name.
#' @param sf a number indicating the scaling factor.
#'        Each coverage value is multiplied by this factor before being reported.
#'        Useful for normalizing coverage by, e.g., reads per million (RPM).
#' @param outfile a character vector indicating the output file name.
#'        If not provided, the result will be internalized in R.
#'
#' @export
#' @return An object of type \code{\link[data.table]{data.table}} with 8 columns,
#'         corresponding to:
#' \enumerate{
#'    \item Chromosome
#'    \item Starting position
#'    \item Ending position
#'    \item Gene name
#'    \item Gene strand (+ or -)
#'    \item Positions. One based positions of the gene.
#'    \item Depth at each position of the gene.
#' }
#'
#' @importFrom utils capture.output
#' @examples
#' data(annotations)
#' bam  <- system.file("extdata", "wt_rep1.bam", package = "yCrypticRNAs")
#' fragments <- bam_to_reads(bam, annotations)
#' coverage(fragments, annotations, sf = 0.069872847)
coverage <- function(reads, annotation, sf = NULL, outfile = NULL){
  #create temp files
  #write bed formatted dataframes to tempfile
  if (is.null(reads))
    stop("The inputfile can't be NULL")
  if(is.null(annotation))
    stop("The annotation can't be NULL")


  if (is.character(reads)){
    a.file = reads
  } else{
    .writeFile(reads, a.file <- tempfile())
  }

  if (is.character(annotation)){
    b.file = annotation
  }else{
    annotation <- as.annotationsSet(annotation)
    .writeFile(annotation, b.file <- tempfile())
  }

  if (is.null(outfile)) {
    file <- tempfile()
  }else{
    file <- outfile
    if (file.info(file)$size != 0)
      unlink(file)
  }

  tryCatch(
    capture.output(
      coverage_cpp(a.file, b.file, c("-d", "-s")),
      file = file, append = T
    ),
    error = function(e) {
      stop(e)
    }
  )

  if (!is.character(reads))
    unlink(a.file)
  if (!is.character(annotation))
    unlink(b.file)

  if(!is.null(sf))
    scale_data(data = file, sf = sf)

  if (!is.null(outfile)) {
    return(paste("see results in file: ", outfile))
  }

  res <- tryCatch(
    data.table::fread(file),
    error = function(e) stop("Error in coverage function: ", e)
  )
  unlink(file)
  return(res)
}
