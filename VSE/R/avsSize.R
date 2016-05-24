#' avsSize
#'
#' This function will output a dataframe of LD size of each tag SNP
#' @param avs A GRanges object which is outputted by loadLd function
#' @keywords AVS,Granges
#' @examples
#' ld<-loadLd(file.path(system.file("extdata", "ld_BCa_raggr.csv", package="VSE")), type="raggr")
#' avs<-makeAVS(ld)
#' avsSize(avs)
#' @importFrom GenomicRanges elementMetadata
#' @export
avsSize <- function(avs){
  no_of_tags <- length(avs)
  avs.ids <- vector()
  size <- vector()
  for (i in 1:no_of_tags){
    avs.ids <- c(avs.ids, as.character(GenomicRanges::elementMetadata(avs[[i]])[1,2]))
    size <- c(size, length(avs[[i]]))
  }
  data.frame(tagID=avs.ids, Size=size)
}
