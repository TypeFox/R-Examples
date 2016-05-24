## $Id: 

if (FALSE){
source("rreport.s")
#source("refManager.s")
#newMarker="mar1";label="label1"; keyword="sae1"
getReference("sae", "ser.adv withdr")
getReference("Osae", "ser.advWITHDR")
getReference("Osae", "ser.advCARDIO")
getReference("sae", "ser.adv cardio")

getRefsByKey("sae")
getLabelsByKey("sae")

getReferenceString("Osae")

getReferenceObject()
}

#' Reference Objects
#'
#' summary
#'
#' details
#'
#' @rdname references
#' @export
#' @examples
#' 1

getReferenceObject <- function(){
  refD <- options("rreport.reference.list")[[1]]
  if (is.null(refD)){
    refD <- data.frame(marker=c(), keyword=c(), label=c())
  }
  refD
}

#' @rdname references
#' @param refD
#' @export
#' @examples
#' 1

putReferenceObject <- function(refD){
  options(rreport.reference.list=refD)
}

#' @rdname references
#' @export
#' @examples
#' 1

print.latexReference <- function(refD){
  if (is.null(refD)){
    cat("The list of markers has not been created. Use function getReferenceObject() to create it\n")
  }else{
    print(refD)
  }
}

#' @rdname references
#' @param newMarker
#' @param keyword
#' @param label
#' @export
#' @examples
#' 1

updateMarkers <- function(newMarker, keyword="", label=""){
  ### puts a new marker into the dataframe of the existing ones
  ### checks if it is different from the existing ones
  ### returns updated latexReference
  refD = getReferenceObject()
  if (newMarker %in% refD$marker){
    stop(paste("Duplicated marker", newMarker))
  }
  newM <-data.frame(marker=newMarker, keyword=keyword, label=label)
  newM$marker <- as.character(newM$marker)
  newM$keyword <- as.character(newM$keyword)
  newM$label <- as.character(newM$label)
  refD = rbind(refD, newM)
  for (n in names(refD)) refD[[n]] <- as.character(refD[[n]])
  putReferenceObject(refD)
}

#' @rdname references
#' @export
#' @examples
#' 1

generateRef <- function(){
  generate <- function(){paste("marker",abs(round(rnorm(1)*(10^8))), sep="")}
  existingMarkers <- getRefsByKey()
  newMarker <- generate()
  while (newMarker %in% existingMarkers){
    newMarker <- generate()
  }
  newMarker
}

#' @rdname references
#' @export
#' @examples
#' 1

getReference <- function(keyword="", label=""){
  newMarker <- generateRef()
  updateMarkers(newMarker = newMarker, keyword=keyword, label=label)
  newMarker
}

#' @rdname references
#' @export
#' @examples
#' 1

getRefsByKey <- function(keyword=NULL){
  ### returns all markers with a given keyword 
  ### if keyword==NULL returns all markers
  refD = getReferenceObject()
  if (!is.null(keyword)){
    refD$marker[refD$keyword==keyword]
  }else{
    refD$marker
  }
}

#' @rdname references
#' @export
#' @examples
#' 1

getLabelsByKey <- function(keyword=NULL){
  ### returns all labels with a given keyword 
  ### if keyword==NULL returns all labels
  refD = getReferenceObject()
  if (!is.null(keyword)){
    refD$label[refD$keyword==keyword]
  }else{
    refD$label
  }
}

#' @rdname references
#' @export
#' @examples
#' 1

getReferenceString <- function(keyword){
  ### returns a vector of strings "see section \\ref{m1} (page\\pageref{m1})"
  ### for all markers with a given keyword 
  markers <- getRefsByKey(keyword)
  labels <- getLabelsByKey(keyword)
  keys <- paste(labels," in section ", "\\ref{", markers, "}", " (page ", "\\pageref{",markers,"}",")", sep="")
  paste("See", paste(keys, collapse=", "))
}

