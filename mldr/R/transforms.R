#' Transformns an MLDR into binary or multiclass datasets
#' @description Transforms an \code{mldr} object into one or serveral binary or multiclass datasets, returning them as \code{data.frame} objects
#' @param mldr The mldr object to transform
#' @param type Indicates the type of transformation to apply. Possible types are:\itemize{
#'  \item \code{"BR"} Produces one or more binary datasets, each one with one label
#'  \item \code{"LP"} Produces a multiclass dataset using each labelset as class label
#'  }
#' @param labels Vector with the label indexes to include in the transformation. All labels will be used if not specified
#' @return A list of data.frames containing the resulting datasets (for BR) or a data.frame with the dataset (for LP).
#' The result is no longer an \code{mldr} object, but a plain \code{data.frame}
#' @examples
#' library(mldr)
#' emotionsbr <- mldr_transform(emotions, type = "BR")
#' emotionslp <- mldr_transform(emotions, type = "LP")
#' @export
mldr_transform <- function(mldr, type = 'BR', labels) {
  if(class(mldr) != 'mldr')
    stop('This method applies only to mldr objects')

  if(missing(labels))
    labels <- mldr$labels$index
  else
    labels <- labels[labels %in% mldr$labels$index]

  if(length(labels) == 0)
    stop('One or more labels have to be selected')

  switch(type,
         BR = mldr_to_BR(mldr, labels = labels),
         LP = mldr_to_LP(mldr, labels = labels)
  )
}

# Internal function to generate BR transformation
mldr_to_BR <- function(mldr, labels) {
  lapply(labels, function(aLabel) {
    binary <- cbind(mldr$dataset[ , mldr$attributesIndexes], mldr$dataset[ , aLabel])
    names(binary)[length(binary)] <- names(mldr$dataset)[aLabel]
    binary
  })
}

# Internal function to generate LP transformation
mldr_to_LP <- function(mldr, labels) {
  cbind(mldr$dataset[,mldr$attributesIndexes],
        classLabel = do.call(paste, c(mldr$dataset[ , labels], sep = "")))
}
