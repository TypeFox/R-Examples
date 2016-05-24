## data_utils.R
##   - Utility functions for R matrices and data frames
##
## RGP - a GP system for R
## 2010-2012 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Patrick Koch, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' Embed columns in a data frame
##'
##' Embeds the columns named \code{cols} in the data frame \code{x} into a space of dimension
##' \code{dimension}.
##'
##' @param x The data frame containing the columns to embed.
##' @param cols A vector a list of the names of the columns to embed.
##' @param dimension The additional dimensions to generate when embedding.
##' @return The data frame, augmented with embedded columns, shortended by \code{dimension} rows.
##'
##' @export
embedDataFrame <- function(x, cols = NULL, dimension = 1) {
  if (dimension == 0)
    return(x)
  d <- dimension + 1
  if (d <= 1 || d >= nrow(x))
    stop("dimension must be greater than 1 and smaller than nrow(x)")
  columns <- if (missing(cols)) names(x) else cols

  embeddings <- lapply(columns, function(column) {
                                  embedded <- as.matrix(embed(x[[column]], d)[,-1])
                                  colnames(embedded) <- paste(column, ".P", 1:ncol(embedded), sep = "")
                                  embedded
                                })
  data.frame(list(x[-1:-(d-1),], embeddings))
}

##' Select a continuous subframe of a data frame
##'
##' Return a continuous subframe of the data frame \code{x} containing \code{size} * \code{nrow(x)}
##' rows from the start, center, or end.
##'
##' @param x The data frame to get a subframe from.
##' @param size The size ratio of the subframe. Must be between 0 and 1.
##' @param pos The position to take the subframe from. Must be \code{"START"}, \code{"CENTER"},
##'   or \code{"END"}.
##' @return A subframe of \code{x}.
##'
##' @export
subDataFrame <- function(x, size = 1.0, pos = "START") {
  if (size < 0 || size > 1) stop("size must be between 0 and 1")
  if (!(pos %in% c("START", "CENTER", "END"))) stop("pos must be either START, CENTER, or END")
  numberOfRows <- nrow(x)
  numberOfSelectedRows <- size * numberOfRows
  if (numberOfSelectedRows == 0) return(NULL)
  startIndex <-
    if ("START" == pos) 1
    else if ("CENTER" == pos) 1 + (numberOfRows - numberOfSelectedRows) %/% 2
    else if ("END" == pos) 1 + numberOfRows - numberOfSelectedRows
    else stop()
  endIndex <- startIndex + numberOfSelectedRows - 1
  x[startIndex:endIndex,]
}
