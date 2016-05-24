#epidase documentation for RData using roxygen2
#put you .RData in the /data directory
#add you roxygen2 description of your data below
#run roxygenize("./pkg") in the root directory of the OutbreakTools package

#' @name FluH1N1pdm2009
#' @title Dataset from the 2009 influenza A/H1N1 pandemic
#' @description This dataset is a \code{list} containing the following objects:
#' \enumerate{
#' \item\code{individuals}: a data frame containing 514 individual IDs as well as their locations.
#' \item\code{dna}: a \code{\link{DNAbin}} object containing 514 genetic sequences of influenza A/H1N1/2009 haemagglutinin (HA).
#' \item\code{tree}: a \code{\link{multiphylo}} object containing the maximum clade credibility (MCC) tree obtained via the Beast analysis of the 514 genetic sequences.
#' }
#' @author Anton Camacho
#' @docType data
#' @references This dataset is part of Trevor Bedford's tutorial on the Beast software: \emph{Inferring spatiotemporal dynamics of the H1N1 influenza pandemic from sequence data},
#' available at \url{https://github.com/trvrb/influenza-dynamics-practical}.
#' @example ../misc/FluH1N1pdmExample.R
#'
NULL
