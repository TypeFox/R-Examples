#' File Converter
#'
#' \code{read4pcadapt} converts \code{.vcf} and \code{.ped} files to an appropriate
#' type of file readable by \code{pcadapt}. You may find the converted file in the
#' current directory.
#'
#' @param input.filename a character string specifying the name of the file to be
#' converted.
#' @param type a character string specifying the type of data to be converted to the
#' \code{pcadapt} format. Supported formats are: \code{ped}, \code{vcf}, \code{lfmm}.
#'
#' @useDynLib pcadapt wrapper_converter
#'
#' @export
#'
read4pcadapt <- function(input.filename,type){
  # input file
  if (class(input.filename) != "character"){
    stop(paste0("File ",input.filename," does not exist."))
  } else if (!file.exists(input.filename)){
    stop(paste0("File ",input.filename," does not exist."))
  }

  if (class(type) != "character" || (!(type %in% c("vcf","ped","lfmm","pcadapt")))){
    stop("Incorrect type.")
  }

  if (type == "ped"){
    .C("wrapper_converter",as.character(input.filename),as.integer(0),PACKAGE = "pcadapt")
  } else if (type == "vcf"){
    .C("wrapper_converter",as.character(input.filename),as.integer(1),PACKAGE = "pcadapt")
  } else if (type == "lfmm"){
    .C("wrapper_converter",as.character(input.filename),as.integer(2),PACKAGE = "pcadapt")
  }
  split.name <- unlist(unlist(strsplit(input.filename, "[.]")))
  
  if ((split.name[-1] %in% c("ped","vcf","lfmm","pcadapt")) && (length(split.name) > 1)){
    aux <- NULL
    for (k in (1:(length(split.name)-1))){
      aux <- paste0(aux,split.name[k])
    }
    aux <- paste0(aux,".pcadapt")
  }
  else {
      aux <- input.filename
  }
  return(aux)
}

#' Population colorization
#'
#' \code{get.score.color} allows the user to display individuals of the same
#' pre-defined population with the same color when using the option
#' \code{"scores"} in \code{pcadapt}.
#'
#' @param pop a list of integers or strings specifying which population the
#' individuals belong to.
#'
#' @examples
#' ## see also ?pcadapt for examples
#'
#' @importFrom grDevices rainbow
#'
#' @keywords internal
#'
#' @export
#'
get.score.color = function(pop){
  pop.split <- list()
  list.ref <- unlist(pop)
  nIND <- length(list.ref)
  idx <- 1
  while (length(list.ref)>0){
    col <- list.ref[1]
    pop.split[[idx]] <- which(pop==col)
    idx <- idx+1
    list.ref <- list.ref[list.ref != col]
  }
  color.list <- rainbow(length(pop.split))
  color.individuals <- array(dim=nIND)
  for (k in 1:length(pop.split)){
    color.individuals[unlist(pop.split[k])] <- color.list[k]
  }
  return(color.individuals)
}


