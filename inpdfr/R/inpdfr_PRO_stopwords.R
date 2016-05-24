#' Load a list of stopwords.
#'
#' \code{getStopWords} returns a list of stopwords.
#'
#' @return A list of vectors with stopwords for French, English, and Spanish languages.
#' @examples
#' getStopWords()
#' @export
getStopWords <- function(){
#   data("exclusionList_FR")
#   data("exclusionList_SP")
#   data("exclusionList_UK")
  exclusionList_FR <- stringi::stri_unescape_unicode(exclusionList_FR)
  exclusionList_SP <- stringi::stri_unescape_unicode(exclusionList_SP)
  Encoding(exclusionList_FR) <- "UTF-8"
  Encoding(exclusionList_SP) <- "UTF-8"
  Encoding(exclusionList_UK) <- "UTF-8"
  return(list(exclusionList_FR, exclusionList_UK, exclusionList_SP))
}

#' Exclude StopWords form the word-occurrence data.frame.
#'
#' Exclude StopWords form the word occurrences data.frame. \code{excludeStopWords}
#'   uses \code{parallel} to perform parallel computation.
#'
#' @param wordF The data.frame containing word occurrences.
#' @param lang The language used ("French", "English", "Spanish").
#' @return The word-occurrence data.frame.
#' @examples
#' \dontrun{
#' excludeStopWords(wordF = myDF, lang = "French")
#' }
#' @export
excludeStopWords <- function(wordF, lang = "English"){
  stopWords <- getStopWords()
  if(lang != "none"){
    if(lang=="French"){
      exclusionList <- stopWords[[1]]
    }else{
      if(lang=="English"){
        exclusionList <- stopWords[[2]]
      }else{
        if(lang=="Spanish"){
          exclusionList <- stopWords[[3]]
        }else{
          if(lang=="user-defined"){
            exclusionList <- stopWords[[4]]
          }
        }
      }
    }
    ncores <- parallel::detectCores()
    if(length(wordF)<ncores){ncores <- length(wordF)}
    cl <- parallel::makeCluster(ncores)
    # parallel::clusterExport(cl = cl, varlist = c("exclusionList"))  ### for testing purposes
    on.exit(parallel::stopCluster(cl))
    part <- parallel::clusterSplit(cl, seq = wordF)
    wordFreqInter <- parallel::parLapply(cl, part, function(partWordFreq){
      for(k in 1:length(partWordFreq)){
        for(i in 1:length(exclusionList)){
          partWordFreq[[k]][[1]] <- partWordFreq[[k]][[1]][partWordFreq[[k]][[1]]$word != exclusionList[i],]
        }
      }
      return(partWordFreq)
    })
    wordF <- unlist(wordFreqInter, recursive = FALSE)
  }
  return(wordF)
}
