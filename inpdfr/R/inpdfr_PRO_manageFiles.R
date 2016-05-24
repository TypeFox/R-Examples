#' List files in a specified directory sorted by extension.
#'
#' List files in a specified directory sorted by extension. The function takes
#'   into account .txt and .pdf files based on \code{strsplit} function.
#'
#' @param mywd A string containing the working directory.
#' @return A list of length 2 with file names sorted by extension (pdf and txt).
#' @examples
#' getListFiles(mywd = getwd())
#' @export
getListFiles <- function(mywd){
  setwd(mywd)
  listFiles <- list.files(pattern = "\\.")
  listFilesExt <- list(pdf = c(), txt = c())
  for(i in listFiles){
    myplitExt <- strsplit(i, split = "\\.")[[1]][2]
    if (myplitExt == "pdf"){listFilesExt$pdf <- c(listFilesExt$pdf, i)}
    if (myplitExt == "txt"){listFilesExt$txt <- c(listFilesExt$txt, i)}
  }
  return(listFilesExt)
}

#' Delete spaces in file names.
#'
#' Delete spaces in file names located in the current working directory.
#'
#' @param vectxt A vector containing character entries corresponding to the names
#'   of files in the current working directory.
#' @return The function returns a logical for each file, with TRUE if the file
#'   has been found, and FALSE otherwise.
#' @examples
#' quitSpaceFromChars(c("my pdf.pdf","my other pdf.pdf"))
#' @export
quitSpaceFromChars <- function(vectxt){ # delete spaces from file names, adapted from : https://gist.github.com/benmarwick/11333467
  isRenamed <- sapply(vectxt, FUN = function(i){
    file.rename(from = i, to = paste0(dirname(i), "/", gsub(" ", "", basename(i))))
  })
  return(isRenamed)
}
